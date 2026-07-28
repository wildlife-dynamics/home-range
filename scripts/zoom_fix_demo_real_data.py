"""Real-data version of zoom_fix_demo.py: fetches Habiba's actual trajectory from
EarthRanger, computes her real ETD home-range percentiles, and renders BEFORE (naive
zoom, assumes a 256px viewport) vs AFTER (the fix in ecoscope's
view_state_from_geodataframes, corrected for a ~1200px viewport) maps side by side.

Reuses bbmm.py's data fetch + render_percentile_map plumbing (render_percentile_map
already calls the real, now-fixed draw_map internally, so it produces the AFTER map
as-is; BEFORE is rendered with an explicit view_state computed via the old formula).
"""

from __future__ import annotations

import logging

import bbmm
import pydeck as pdk

from ecoscope.platform.tasks.config import EtdMethodArgs, call_home_range_from_args
from ecoscope.platform.tasks.preprocessing import convert_trajectory_to_relocations
from ecoscope.platform.tasks.results import create_geojson_layer, draw_map
from ecoscope.platform.tasks.results._map_utils import TileLayer
from ecoscope.platform.tasks.results._pydeck import GeoJSONLayerStyle, LegendFromDataframe, LegendStyle
from ecoscope.platform.tasks.transformation import apply_color_map, convert_column_values_to_string, convert_crs, map_columns

logger = logging.getLogger(__name__)

SUBJECT_NAME = "Habiba"
SUBJECT_ID = bbmm.SUBJECTS[SUBJECT_NAME]


def render_with_explicit_view_state(percentiles_gdf, view_state, output_path: str, title: str) -> None:
    """Same rendering as bbmm.render_percentile_map, but with an explicit view_state
    override (so we can force the OLD, unfixed naive zoom for the BEFORE comparison)."""
    df = convert_column_values_to_string(percentiles_gdf, columns=["percentile"])
    df = apply_color_map(df, input_column_name="percentile", colormap="RdYlGn_r", output_column_name="percentile_color")
    df = convert_crs(df, crs="EPSG:4326")
    df = map_columns(
        df,
        drop_columns=[],
        retain_columns=[],
        rename_columns={"percentile": "Percentile", "area_sqkm": "Area"},
        raise_if_not_found=True,
    )

    layer = create_geojson_layer(
        geodataframe=df,
        layer_style=GeoJSONLayerStyle(
            filled=True,
            stroked=False,
            get_fill_color="percentile_color",
            opacity=0.7,
        ),
        legend=LegendFromDataframe(
            title="Home Range %",
            label_column="Percentile",
            label_suffix=" %",
            color_column="percentile_color",
            sort="ascending",
        ),
        tooltip_columns=["Percentile", "Area"],
        zoom=False,
    )
    tile_layer = TileLayer(layer_name="TERRAIN", opacity=1.0)
    html = draw_map(
        geo_layers=[layer],
        tile_layers=[tile_layer],
        static=False,
        output_type="html",
        view_state=view_state,
        max_zoom=15,
        title=title,
        legend_style=LegendStyle(placement="bottom-right"),
    )
    with open(output_path, "w") as f:
        f.write(html)
    logger.info(f"wrote {output_path}  ->  {view_state}")


def main() -> None:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    traj_gdf = bbmm.fetch_subject_trajectory(client, SUBJECT_ID)
    logger.info(f"{SUBJECT_NAME}: {len(traj_gdf)} trajectory segments fetched")

    args = EtdMethodArgs(
        opacity=0.7,
        crs="EPSG:3857",
        max_speed_factor=1.05,
        expansion_factor=1.3,
        percentiles=[50.0, 60.0, 70.0, 80.0, 90.0, 99.999],
    )
    result = call_home_range_from_args(traj_gdf, args)
    logger.info(f"{SUBJECT_NAME}: computed {len(result)} ETD percentile polygons")

    # AFTER: render_percentile_map calls draw_map with view_state=None internally,
    # which now goes through the FIXED view_state_from_geodataframes.
    bbmm.render_percentile_map(
        result,
        "/tmp/habiba_real_zoomfix_after.html",
        colorize=True,
        filled=True,
        stroked=False,
        opacity=0.7,
        legend_title="Home Range %",
    )

    # BEFORE: reproduce the OLD (unfixed) naive zoom explicitly for comparison.
    wgs84 = result.to_crs("EPSG:4326")
    bounds = wgs84.total_bounds
    bbox = [[bounds[0], bounds[1]], [bounds[2], bounds[3]]]
    naive_zoom = pdk.data_utils.viewport_helpers.bbox_to_zoom_level(bbox)
    center_lon = (bounds[0] + bounds[2]) / 2
    center_lat = (bounds[1] + bounds[3]) / 2
    before_view_state = pdk.ViewState(longitude=center_lon, latitude=center_lat, zoom=naive_zoom)

    render_with_explicit_view_state(
        result,
        before_view_state,
        "/tmp/habiba_real_zoomfix_before.html",
        "BEFORE fix (naive, zoomed out)",
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    main()
