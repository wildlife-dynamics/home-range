"""Classic Brownian Bridge Movement Model (BBMM) home-range estimation.

Standalone prototype - not wired into spec.yaml/wt-compiler. Implements the
original (non-dynamic) BBMM from Horne, Garton, Krone & Lewis (2007),
"Analyzing animal movements using Brownian bridges", Ecology 88(9):2354-63.
The bridge-variance and log-likelihood formulas match both the reference R
implementation (CRAN package `BBMM`) and Eqn 1 of Kranstauber et al. (2012,
J. Animal Ecology 81:738-746), which reproduces Horne et al.'s original
likelihood verbatim as the baseline for their dynamic BBMM extension.
Locations are required to be in an equal-area projection (per Kranstauber
et al. 2012, Table 1); the motion-variance MLE uses only every *other*
interior fix, so the leave-one-out tests are statistically independent.

Model summary
-------------
Between each pair of consecutive fixes (z_i, z_i+1) separated by time_lag,
the animal's true position at elapsed fraction alpha in [0, 1] is modeled as
bivariate normal:

    mean(alpha)     = (1 - alpha) * z_i + alpha * z_i+1
    variance(alpha) = time_lag * alpha * (1 - alpha) * sigma_m2
                       + (1 - alpha)^2 * location_error^2
                       + alpha^2 * location_error^2

sigma_m2 (the Brownian motion variance, i.e. the animal's mobility) is a
single value for the whole track, estimated by maximum likelihood: treat
each interior fix as if it were unobserved and ask how likely its actual
position is under the bridge predicted by its two neighbours.

The utilization distribution (UD) is the sum, over every segment, of that
segment's bridge density integrated over elapsed time - approximated here
by summing the density at discrete time steps and scaling by the segment's
real duration.

Reuses ecoscope's own raster/percentile conventions (ecoscope.io.raster,
ecoscope.analysis.percentile.get_percentile_area) so output shapes match
the ETD home-range pipeline's GeoTIFF + percentile-table pair.
"""

from __future__ import annotations

import logging
import os

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from ecoscope.analysis.percentile import get_percentile_area
from ecoscope.io import raster
from ecoscope.platform.connections import EarthRangerConnection
from ecoscope.relocations import Relocations
from ecoscope.trajectory import Trajectory

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ── Configuration ────────────────────────────────────────────────────────────

DATA_SOURCE = "mep"
SUBJECT_GROUP = "Elephants"
SUBJECTS = {
    "Habiba": "64444ed7-72ec-4531-a2b1-fb25c7197b2d",
    "Salif Keita": "b8be28f7-8c20-46d9-85a5-fd817351bde5",
}

CRS = "ESRI:102022"  # Africa Albers Equal Area Conic - Horne et al. 2007 requires an equal-area projection
LOCATION_ERROR_METERS = 20.0  # typical GPS collar accuracy
TIME_STEP_SECONDS = 60.0  # bridge-integration discretization interval
MAX_STEPS_PER_SEGMENT = 50  # cap per-segment steps regardless of gap size (large gaps -> very diffuse bridge anyway)
GRID_SCALE_FACTOR = 500  # matches ecoscope.analysis.UD.grid_size_from_geographic_extent's own default usage in ETD
EXPANSION_FACTOR = 1.3  # matches ETD's own default grid padding
PERCENTILES = [50.0, 60.0, 70.0, 80.0, 90.0, 95.0, 99.999]
MAX_BRIDGE_RADIUS_SIGMA = 4.0  # restrict density computation to +/- N sigma around each bridge, for speed
MAX_OUTLIER_DISTANCE_KM = 300.0  # drop fixes farther than this from the median position (GPS error rejection)

OUTPUT_DIR = os.environ.get("BBMM_OUTPUT_DIR", "/tmp/workflows/bbmm/output")


# ── Data fetch ───────────────────────────────────────────────────────────────


def _drop_gps_outliers(gdf: gpd.GeoDataFrame, max_distance_km: float = MAX_OUTLIER_DISTANCE_KM) -> gpd.GeoDataFrame:
    """Drop fixes implausibly far (in degrees-as-proxy for km) from the median position.

    Neither `filter="clean"` nor the standard junk-coordinate filter ([180,90]/[0,0]/[1,1])
    catches isolated GPS-error fixes at arbitrary bad coordinates - e.g. Habiba's real data
    includes two fixes ~4000km from the rest of her cluster with no special-case coordinate
    value. A simple distance-from-median-position cutoff catches this class of error directly.
    """
    median_lon, median_lat = gdf.geometry.x.median(), gdf.geometry.y.median()
    # ~111 km per degree of latitude; close enough for a coarse outlier cutoff
    dist_km = np.sqrt((gdf.geometry.x - median_lon) ** 2 + (gdf.geometry.y - median_lat) ** 2) * 111.0
    keep = dist_km <= max_distance_km
    n_dropped = (~keep).sum()
    if n_dropped:
        logger.warning(f"Dropping {n_dropped} GPS outlier fix(es) more than {max_distance_km} km from the median position")
    return gdf[keep].copy()


def fetch_subject_trajectory(client, subject_id: str) -> gpd.GeoDataFrame:
    """Fetch and build a single subject's trajectory (segment_start/segment_end etc.) from raw relocations."""
    relocs = client.get_subjectgroup_observations(subject_group_name=SUBJECT_GROUP, filter="clean")
    gdf = relocs.gdf
    subject_gdf = gdf[gdf["groupby_col"] == subject_id].copy()
    if subject_gdf.empty:
        raise ValueError(f"No observations found for subject_id={subject_id!r}")
    subject_gdf = _drop_gps_outliers(subject_gdf)

    traj = Trajectory.from_relocations(Relocations(gdf=subject_gdf))
    traj_gdf = traj.gdf.to_crs(CRS)
    return traj_gdf.sort_values("segment_start").reset_index(drop=True)


def relocations_from_trajectory(traj_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Reconstruct a point-relocations GeoDataFrame (geometry + fixtime) from a
    Trajectory's LineString segments, via _extract_points. This is lossless - n
    fixes become n-1 segments and back again exactly, per Trajectory.from_relocations'
    own _create_multitraj (each segment is literally one pair of consecutive fixes) -
    so this recovers the same n original fixes, just without any segment-derived
    columns (segment_start/segment_end/etc.) still attached.

    Useful for methods like MCP that don't need segment/time-lag structure at all
    (see MCPRange.cs: CalculateMCPRange takes a Trajectory parameter, but only ever
    reads thePath.DataPoints.RelocationsArray - the raw relocation list) but still
    want to reuse fetch_subject_trajectory's single fetch/filter/reproject path.
    """
    xy, t = _extract_points(traj_gdf)
    geometry = gpd.points_from_xy(xy[:, 0], xy[:, 1])
    fixtime = pd.to_datetime(t, unit="s", utc=True)
    return gpd.GeoDataFrame({"fixtime": fixtime, "geometry": geometry}, crs=traj_gdf.crs)


# ── Brownian motion variance (MLE) ──────────────────────────────────────────


def _bridge_variance(time_lag: float, alpha: np.ndarray, sigma_m2: float, location_error: float) -> np.ndarray:
    return time_lag * alpha * (1 - alpha) * sigma_m2 + ((1 - alpha) ** 2 + alpha**2) * location_error**2


def _extract_points(traj_gdf: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Reconstruct the point sequence (xy, t) from a trajectory's LineString segments.

    Each row's geometry is a 2-point LineString from segment_start to segment_end;
    consecutive segments share an endpoint, so the full point path is the first
    segment's start point followed by every segment's end point.
    """
    starts = np.array([geom.coords[0] for geom in traj_gdf.geometry])
    ends = np.array([geom.coords[-1] for geom in traj_gdf.geometry])
    xy = np.vstack([starts[[0]], ends])
    t = (
        pd.concat([traj_gdf["segment_start"].iloc[[0]], traj_gdf["segment_end"]])
        .astype("int64")
        .to_numpy()
        / 1e9
    )
    return xy, t


def _neg_log_likelihood(sigma_m2: float, midpoints: pd.DataFrame, location_error: float) -> float:
    variance = _bridge_variance(midpoints["jump_time"].to_numpy(), midpoints["alpha"].to_numpy(), sigma_m2, location_error)
    d2 = midpoints["d2"].to_numpy()
    log_lik = -np.log(2 * np.pi * variance) - d2 / (2 * variance)
    return -log_lik.sum()


def compute_midpoints(traj_gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    """Odd-i leave-one-out test points (Horne et al. 2007, Eqn 1), each tagged with the
    two segment indices (into the full per-segment array used by compute_bbmm_ud) it
    borders. Only every *other* interior fix is used, so the leave-one-out tests are
    statistically independent - consecutive interior points share an overlapping
    neighbour, which would violate the i.i.d. assumption behind the likelihood product.
    Shared by both the classic (single sigma_m^2) and dynamic (windowed) MLE.
    """
    xy, t = _extract_points(traj_gdf)

    z_prev, z_curr, z_next = xy[:-2], xy[1:-1], xy[2:]
    t_prev, t_curr, t_next = t[:-2], t[1:-1], t[2:]
    segment_before = np.arange(len(xy) - 2)  # segment index just before z_curr
    segment_after = segment_before + 1  # segment index just after z_curr

    # Subsample to only odd i (position 0 in these arrays == i=1, so [0::2] selects i=1,3,5,...)
    sl = slice(0, None, 2)
    z_prev, z_curr, z_next = z_prev[sl], z_curr[sl], z_next[sl]
    t_prev, t_curr, t_next = t_prev[sl], t_curr[sl], t_next[sl]
    segment_before, segment_after = segment_before[sl], segment_after[sl]

    jump_time = t_next - t_prev
    valid = jump_time > 0
    alpha = np.where(valid, (t_curr - t_prev) / np.where(valid, jump_time, 1.0), np.nan)
    mu = (1 - alpha[:, None]) * z_prev + alpha[:, None] * z_next
    d2 = np.sum((z_curr - mu) ** 2, axis=1)

    midpoints = pd.DataFrame(
        {
            "jump_time": jump_time,
            "alpha": alpha,
            "d2": d2,
            "segment_before": segment_before,
            "segment_after": segment_after,
        }
    )
    return midpoints[valid & midpoints["alpha"].between(0, 1)].reset_index(drop=True)


def estimate_motion_variance(traj_gdf: gpd.GeoDataFrame, location_error: float) -> float:
    """MLE for a single, whole-track sigma_m^2 (Horne et al. 2007, Eqn 1)."""
    midpoints = compute_midpoints(traj_gdf)
    if midpoints.empty:
        raise ValueError("Not enough interior fixes with valid time gaps to estimate motion variance.")

    result = minimize_scalar(
        _neg_log_likelihood,
        bounds=(1.0, 1_000_000.0),
        method="bounded",
        args=(midpoints, location_error),
    )
    logger.info(f"Estimated Brownian motion variance (sigma_m^2) = {result.x:.2f} m^2/s (from {len(midpoints)} fixes)")
    return float(result.x)


# ── UD raster computation ───────────────────────────────────────────────────


def _build_grid(traj_gdf: gpd.GeoDataFrame, pixel_size: float, expansion_factor: float):
    x_min, y_min, x_max, y_max = traj_gdf.geometry.total_bounds
    if expansion_factor > 1.0:
        dx = (x_max - x_min) * (expansion_factor - 1.0) / 2.0
        dy = (y_max - y_min) * (expansion_factor - 1.0) / 2.0
        x_min, x_max = x_min - dx, x_max + dx
        y_min, y_max = y_min - dy, y_max + dy

    raster_profile = raster.RasterProfile(
        pixel_size=pixel_size,
        crs=CRS,
        nodata_value=np.nan,
        band_count=1,
        raster_extent=raster.RasterExtent(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max),
    )
    num_rows, num_cols = raster_profile.rows, raster_profile.columns

    col_centers = x_min + pixel_size * (np.arange(num_cols) + 0.5)
    row_centers = y_max - pixel_size * (np.arange(num_rows) + 0.5)
    return raster_profile, col_centers, row_centers


def compute_bbmm_ud(
    traj_gdf: gpd.GeoDataFrame,
    sigma_m2: float | np.ndarray,
    location_error: float,
    pixel_size: float,
    time_step_seconds: float = TIME_STEP_SECONDS,
    expansion_factor: float = EXPANSION_FACTOR,
    max_steps_per_segment: int = MAX_STEPS_PER_SEGMENT,
    window_padding_sigma: float = MAX_BRIDGE_RADIUS_SIGMA,
    max_data_gap_seconds: float | None = None,
) -> raster.RasterData:
    """Compute the BBMM UD raster. `sigma_m2` may be a single value (classic BBMM) or a
    per-segment array of length len(traj_gdf) - 1 (dynamic BBMM, one value per segment).

    `max_steps_per_segment` caps how finely each segment's bridge-integration is
    discretized, regardless of its own time gap (a large gap otherwise produces a
    very diffuse bridge anyway, so extra steps buy little). `window_padding_sigma`
    sizes the padded rectangular window each segment's density is computed
    over (a vectorization/performance device, not part of the Horne et al. 2007
    formula itself) - both default to this module's own constants, unchanged.

    `max_data_gap_seconds` (default None: no exclusion) drops segments with a time
    gap >= this threshold entirely, rather than modeling them with a Brownian
    bridge - bridge variance scales linearly with time_lag, so one atypically long
    gap (e.g. a multi-day collar dropout amid otherwise hourly fixes) produces a
    hugely diffuse, near-uniform blob dominating the whole UD, even though it's
    really just "we have no data here", not meaningful movement uncertainty.
    ArcMET's own BBMM implementation defaults this to 14400 (4 hours)."""
    raster_profile, col_centers, row_centers = _build_grid(traj_gdf, pixel_size, expansion_factor)
    num_rows, num_cols = raster_profile.rows, raster_profile.columns
    ud = np.zeros((num_rows, num_cols), dtype=np.float64)

    xy, t = _extract_points(traj_gdf)

    n_segments = len(xy) - 1
    sigma_m2_per_segment = np.broadcast_to(np.asarray(sigma_m2, dtype=np.float64), (n_segments,))

    for i in range(n_segments):
        z0, z1 = xy[i], xy[i + 1]
        time_lag = t[i + 1] - t[i]
        if time_lag <= 0:
            continue
        if max_data_gap_seconds is not None and time_lag >= max_data_gap_seconds:
            continue

        n_steps = min(max_steps_per_segment, max(2, int(np.ceil(time_lag / time_step_seconds)) + 1))
        alphas = np.linspace(0.0, 1.0, n_steps)
        dt = time_lag / (n_steps - 1)

        variances = _bridge_variance(time_lag, alphas, sigma_m2_per_segment[i], location_error)
        max_sigma = np.sqrt(variances.max())
        pad = window_padding_sigma * max_sigma + pixel_size

        seg_x_min, seg_x_max = min(z0[0], z1[0]) - pad, max(z0[0], z1[0]) + pad
        seg_y_min, seg_y_max = min(z0[1], z1[1]) - pad, max(z0[1], z1[1]) + pad

        col_mask = (col_centers >= seg_x_min) & (col_centers <= seg_x_max)
        row_mask = (row_centers >= seg_y_min) & (row_centers <= seg_y_max)
        if not col_mask.any() or not row_mask.any():
            continue

        local_x = col_centers[col_mask]
        local_y = row_centers[row_mask]
        X, Y = np.meshgrid(local_x, local_y)

        local_accum = np.zeros_like(X)
        for alpha, variance in zip(alphas, variances):
            mu_x = (1 - alpha) * z0[0] + alpha * z1[0]
            mu_y = (1 - alpha) * z0[1] + alpha * z1[1]
            ztz = (X - mu_x) ** 2 + (Y - mu_y) ** 2
            local_accum += (1.0 / (2 * np.pi * variance)) * np.exp(-ztz / (2 * variance))

        local_accum *= dt  # Riemann-sum approximation of the time integral over this segment
        row_idx = np.where(row_mask)[0]
        col_idx = np.where(col_mask)[0]
        ud[np.ix_(row_idx, col_idx)] += local_accum

        if (i + 1) % 5000 == 0:
            logger.info(f"  ...processed {i + 1}/{n_segments} segments")

    cell_area = pixel_size * pixel_size
    total_mass = ud.sum() * cell_area
    if total_mass > 0:
        ud /= total_mass

    return raster.RasterData(data=ud.astype("float32"), crs=raster_profile.crs, transform=raster_profile.transform), raster_profile


# ── Map rendering (mirrors the ETD pipeline's exact recipe) ─────────────────


def _to_nonoverlapping_rings(df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Convert nested/cumulative percentile polygons - each larger percentile's
    polygon fully contains all smaller ones' (convex hull is monotonic under set
    inclusion, and MCP's percentile point-subsets are themselves nested) - into
    non-overlapping annular rings for display: each ring is that percentile's
    polygon minus the next-smaller percentile's polygon.

    Needed because stacking N overlapping semi-transparent layers compounds their
    visible opacity to 1-(1-opacity)^N (e.g. 7 nested layers at opacity=0.25 stack
    to ~0.87, nearly solid) even though each individual layer's opacity looks low.
    Non-overlapping rings render at exactly the requested opacity everywhere.
    Properties (e.g. the cumulative MCPArea/"Area" tooltip value) are left as-is -
    only the display geometry changes, not what it's labeled with.

    Differencing must walk smallest-to-largest, but apply_color_map's categorical
    branch assigns colors by row *encounter order* (not by value), and the rest of
    this pipeline expects descending order (largest first) for that to come out
    red=smallest/green=largest - so the ascending working order is reversed back
    to descending before returning.
    """
    ascending = df.sort_values("percentile").reset_index(drop=True).copy()
    prev_geom = None
    rings = []
    for geom in ascending.geometry:
        rings.append(geom.difference(prev_geom) if prev_geom is not None else geom)
        prev_geom = geom
    ascending["geometry"] = rings
    return ascending.iloc[::-1].reset_index(drop=True)


def render_percentile_map(
    percentiles_gdf: gpd.GeoDataFrame,
    output_path: str,
    colorize: bool = True,
    filled: bool = True,
    stroked: bool = True,
    opacity: float = 0.6,
    rings: bool = False,
    track_gdf: gpd.GeoDataFrame | None = None,
    track_color: list[int] = [0, 102, 204],
    track_opacity: float = 0.7,
    track_point_radius: float = 3,
    track_on_top: bool = True,
    legend_title: str = "Time Spent",
) -> None:
    """colorize=False renders unfilled outlines only (no RdYlGn_r fill/legend) -
    useful to inspect polygon shapes/boundaries without the percentile color coding.
    When colorize=True, each polygon is (if filled=True) filled with its RdYlGn_r
    percentile color at the given opacity (an rgba fill from apply_color_map), and
    (if stroked=True) outlined in black - stroked=False matches the ETD production
    pipeline's own etd_pct_layer config (spec.yaml), which has no outline at all.
    filled=False renders a completely transparent (unfilled) polygon - just a
    percentile-colored stroke - useful with track_gdf so the track is fully visible
    with nothing drawn over it; the legend still reflects the same percentile_color
    column either way.

    rings=True converts nested/cumulative percentile polygons into non-overlapping
    display rings first (see _to_nonoverlapping_rings) - use this whenever the
    polygons are cumulative (each larger percentile contains the smaller ones, as
    MCP's do) and opacity should look uniform rather than compounding in the core.

    track_gdf, if given, overlays its geometries (raw relocation Points, or a
    trajectory's LineString segments) in track_color at track_opacity, with
    track_point_radius (pixels) for Points - Points render as small dots,
    LineStrings as a connecting line - so the actual GPS fixes/movement path are
    visible alongside the home-range polygons.
    track_on_top controls z-order: True (default) draws it above the percentile
    layer, False draws it underneath (so the polygon fill/stroke draws over it -
    pair with a low track_opacity to still see it through the fill).
    """
    from ecoscope.platform.tasks.results import create_geojson_layer, draw_map
    from ecoscope.platform.tasks.results._map_utils import TileLayer
    from ecoscope.platform.tasks.results._pydeck import GeoJSONLayerStyle, LegendFromDataframe, LegendStyle
    from ecoscope.platform.tasks.transformation import apply_color_map, convert_column_values_to_string, convert_crs, map_columns

    df = percentiles_gdf
    if rings:
        df = _to_nonoverlapping_rings(df)
    df = convert_column_values_to_string(df, columns=["percentile"])
    if colorize:
        df = apply_color_map(df, input_column_name="percentile", colormap="RdYlGn_r", output_column_name="percentile_color")
    df = convert_crs(df, crs="EPSG:4326")
    df = map_columns(
        df,
        drop_columns=[],
        retain_columns=[],
        rename_columns={"percentile": "Percentile", "area_sqkm": "Area"},
        raise_if_not_found=True,
    )

    if colorize:
        layer_style = GeoJSONLayerStyle(
            filled=filled,
            stroked=stroked,
            get_fill_color="percentile_color" if filled else None,
            get_line_color=[0, 0, 0] if filled else "percentile_color",
            get_line_width=2 if stroked else 0,
            line_width_min_pixels=1 if stroked else 0,
            opacity=opacity,
        )
        legend = LegendFromDataframe(
            title=legend_title,
            label_column="Percentile",
            label_suffix=" %",
            color_column="percentile_color",
            sort="ascending",
        )
    else:
        layer_style = GeoJSONLayerStyle(
            filled=False,
            stroked=True,
            get_line_color=[0, 0, 0],
            get_line_width=2,
        )
        legend = None

    layer = create_geojson_layer(
        geodataframe=df,
        layer_style=layer_style,
        legend=legend,
        tooltip_columns=["Percentile", "Area"],
        zoom=True,
    )
    # track_on_top controls draw order: whichever layer is appended last renders on top.
    geo_layers = []
    track_layer = None
    if track_gdf is not None:
        track_df = convert_crs(track_gdf[["geometry"]], crs="EPSG:4326")
        track_layer = create_geojson_layer(
            geodataframe=track_df,
            layer_style=GeoJSONLayerStyle(
                filled=True,
                stroked=True,
                get_fill_color=track_color,
                get_line_color=track_color,
                get_line_width=1,
                line_width_min_pixels=1,
                get_point_radius=track_point_radius,
                point_radius_units="pixels",
                point_radius_min_pixels=1,
                opacity=track_opacity,
            ),
            zoom=False,
        )
        if not track_on_top:
            geo_layers.append(track_layer)
    geo_layers.append(layer)
    if track_gdf is not None and track_on_top:
        geo_layers.append(track_layer)
    tile_layer = TileLayer(
        # TileLayer's own model_validator overwrites `url` whenever `layer_name` matches
        # a known preset (see ecoscope.platform.tasks.results._map_utils), so a custom
        # url= is silently ignored here - "OpenStreetMap" forces tile.openstreetmap.org,
        # whose usage policy blocks unidentified automated/app traffic (403). Use the
        # "TERRAIN" preset instead, which already resolves to the ArcGIS World_Topo_Map
        # tile source used elsewhere in this repo.
        layer_name="TERRAIN",
        opacity=1.0,
    )
    html = draw_map(
        geo_layers=geo_layers,
        tile_layers=[tile_layer],
        static=False,
        output_type="html",
        view_state=None,
        max_zoom=15,
        legend_style=LegendStyle(placement="bottom-right"),
    )
    with open(output_path, "w") as f:
        f.write(html)
    logger.info(f"Wrote map: {output_path}")


# ── Main ─────────────────────────────────────────────────────────────────────


def run_subject(client, subject_name: str, subject_id: str, output_dir: str) -> None:
    logger.info(f"=== {subject_name} ===")
    traj_gdf = fetch_subject_trajectory(client, subject_id)
    logger.info(f"{len(traj_gdf)} trajectory segments fetched")

    sigma_m2 = estimate_motion_variance(traj_gdf, LOCATION_ERROR_METERS)

    from ecoscope.analysis.UD import grid_size_from_geographic_extent

    pixel_size = grid_size_from_geographic_extent(traj_gdf, scale_factor=GRID_SCALE_FACTOR)
    logger.info(f"Grid pixel size: {pixel_size} m")

    raster_data, raster_profile = compute_bbmm_ud(traj_gdf, sigma_m2, LOCATION_ERROR_METERS, pixel_size)

    slug = subject_name.lower().replace(" ", "_")
    os.makedirs(output_dir, exist_ok=True)
    raster_path = os.path.join(output_dir, f"{slug}_bbmm_raster.tif")

    write_array = raster_data.data.copy()
    write_array[np.isnan(write_array) | (write_array == 0)] = raster_profile.nodata_value
    raster.RasterPy.write(write_array, fp=raster_path, **raster_profile)
    logger.info(f"Wrote GeoTIFF: {raster_path}")

    percentiles_gdf = get_percentile_area(PERCENTILES, raster_data, subject_id=subject_name)
    percentiles_gdf = percentiles_gdf.set_geometry("geometry", crs=raster_profile.crs)
    percentiles_gdf["area_sqkm"] = percentiles_gdf.geometry.area / 1e6

    csv_path = os.path.join(output_dir, f"{slug}_bbmm_percentiles.csv")
    percentiles_gdf.drop(columns="geometry").to_csv(csv_path, index=False)
    logger.info(f"Wrote percentile table: {csv_path}")

    for _, row in percentiles_gdf.sort_values("percentile").iterrows():
        logger.info(f"  {row['percentile']:>7.3f}% isopleth area: {row['area_sqkm']:.2f} km^2")

    map_path = os.path.join(output_dir, f"{slug}_bbmm_map.html")
    # stroked=False: BBMM is a density-surface method like ETD, not a hull method
    # like MCP - many small isopleth fragments look like a messy quilted mesh with
    # an outline on every one (matches ETD's own production etd_pct_layer config).
    render_percentile_map(percentiles_gdf, map_path, stroked=False)


def main() -> None:
    client = EarthRangerConnection.client_from_named_connection(DATA_SOURCE)
    for subject_name, subject_id in SUBJECTS.items():
        run_subject(client, subject_name, subject_id, OUTPUT_DIR)


if __name__ == "__main__":
    main()
