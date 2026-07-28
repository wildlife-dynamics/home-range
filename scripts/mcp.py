"""Minimum Convex Polygon (MCP) home-range estimation.

A line-for-line port of ArcMET's own C# implementation at
/Users/zak/Documents/w-dynamics/arcmet/ArcMET_Base/ArcMET_Base/MCPRange.cs
(class MCPRange.CalculateMCPRange), not just the same textbook method -
every numeric step below is traceable to a specific line in that file:

  1. Centroid: add every fix to an ESRI Multipoint and take its CentroidEx
     (MCPRange.cs:259-267) - the plain arithmetic mean of all fix
     coordinates, duplicates included (no dedup).
  2. Distance: Calculations.SeparationBetweenPointsMeters (MCPRange.cs:275).
     For a *projected* spatial reference that function (Calculations.cs:155)
     takes the straight-line (Euclidean) distance between the two points and
     multiplies by the CRS's MetersPerUnit - i.e. plain Euclidean distance
     in projected meters when, as here, the CRS's linear unit already is
     meters. (Its geographic-CRS branch, an inverse-geodesic calculation, is
     never hit because bbmm.fetch_subject_trajectory always reprojects to
     bbmm.CRS, a projected equal-area CRS in meters.)
  3. Sort fixes by that distance, ascending (MCPRange.cs:282).
  4. cutoff = floor(p * n) fixes are kept, the closest ones first
     (MCPRange.cs:287, 290-294); percentiles are themselves iterated in
     ascending order, matching MCPRangeForm.cs:109's `percentiles.Sort()`.
  5. actualPercentile = 100 * floor(p*n) / n (MCPRange.cs:286).
  6. Convex hull of the retained fixes via ITopologicalOperator2.ConvexHull
     (MCPRange.cs:298-299) == shapely's MultiPoint(...).convex_hull.
  7. MCPArea = IArea.Area on that hull (MCPRange.cs:302-303) - raw area in
     whatever units the CRS uses, m^2 here since bbmm.CRS is meters-based.
     No km^2 conversion: ArcMET never does one either.
  8. Per-percentile calculation errors (e.g. fewer than 3 points at that
     percentile, MCPRange.cs:253-254) are caught and recorded rather than
     aborting the whole run (MCPRange.cs:344-348), matching ArcMET's
     Errors-list-per-calculation behaviour; other percentiles still compute.

Output schema matches MCPRange.cs's DefineFields() (MCPRange.cs:396-474)
exactly - MovDataID, CalcID, StartDate, EndDate, ChosenPercentile,
ActualPercentile, MCPArea, in that column order - plus a second table
mirroring MCPRangeMethodResults.WriteCalcSummaryTable's summary fields
(MCPRange.cs:34-98): MovDataID, CalcID, CalcInitiated, CalcFinished,
StartDate, EndDate, ResultFCName, Errors.

Unlike BBMM/dBBMM/aKDE, MCP produces no continuous utilization-distribution
surface (no GeoTIFF) - it's a purely geometric method with no probabilistic
movement model behind it, so there's nothing to rasterize. It's also the
oldest and now most heavily criticized home-range estimator (sensitive to
sample size and outliers, and it includes area the animal demonstrably
never visited, e.g. across a lake) - implemented here for comparison
against BBMM/dBBMM/aKDE, not as a recommended first choice.

Fetches via bbmm.fetch_subject_trajectory (the same single fetch/filter/reproject
path BBMM/dBBMM/aKDE use), then reconstructs the raw point relocations from that
Trajectory's segments via bbmm.relocations_from_trajectory. MCP's own algorithm
never needs the segment/time-lag structure itself - only the fix positions (plus
their timestamps' min/max for StartDate/EndDate) - matching ArcMET's own
MCPRange.cs: CalculateMCPRange takes a Trajectory parameter, but only ever reads
thePath.DataPoints.RelocationsArray, the raw relocation list, never anything
segment-derived.

Reuses bbmm.py's data fetch, CRS, and map plumbing.
"""

from __future__ import annotations

import datetime as dt
import logging
import os

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import MultiPoint

import bbmm

logger = logging.getLogger(__name__)

# ── Configuration ────────────────────────────────────────────────────────────
# ChosenPercentile values as fractions of 1 (MCPRangeForm.cs:104: txtPercentage
# is parsed as "value / 100.0"). Iterated in ascending order to match
# MCPRangeForm.cs:109's percentiles.Sort() (ascending, no reverse).
# PERCENTILES = sorted(p / 100.0 for p in bbmm.PERCENTILES)
# Trial run: just 50/80/99 instead of the full bbmm.PERCENTILES set.
PERCENTILES = sorted(p / 100.0 for p in [50.0, 80.0, 99.0])

CALC_ID = "MCP"  # ArcMET's CalcID is a free-text field the user types into the form (txtCalcID); fixed here for a scripted run.

OUTPUT_DIR = os.environ.get("MCP_OUTPUT_DIR", "/tmp/workflows/mcp/output")


def compute_mcp(relocs_gdf: gpd.GeoDataFrame, percentiles: list[float] = PERCENTILES) -> tuple[gpd.GeoDataFrame, list[str]]:
    """Percentile MCP, replicating MCPRange.CalculateMCPRange (MCPRange.cs:195-383).

    `relocs_gdf` is a plain point GeoDataFrame (bbmm.fetch_subject_relocations) -
    order doesn't matter to the algorithm itself (centroid/distance/hull are all
    order-independent), matching ArcMET's own MCPRange.cs which never relies on
    fix ordering either.

    `percentiles` are fractions of 1 (e.g. 0.95 for 95%), matching ArcMET's
    own internal representation after MCPRangeForm.cs:104 divides the user's
    percent input by 100.

    Returns (rows_gdf, errors) - errors is a list of per-percentile failure
    messages, mirroring MCPRangeMethodResults.Errors (MCPRange.cs:344-348):
    a failure on one percentile is recorded but does not stop the others.
    """
    xy = np.column_stack([relocs_gdf.geometry.x.to_numpy(), relocs_gdf.geometry.y.to_numpy()])
    n = len(xy)

    # MCPRange.cs:259-267: Multipoint.CentroidEx == arithmetic mean of every
    # point added to the multipoint (all n fixes, duplicates included).
    centroid = xy.mean(axis=0)

    # MCPRange.cs:273-279: Calculations.SeparationBetweenPointsMeters per fix.
    # For a projected CRS in meters this is plain Euclidean distance
    # (Calculations.cs:155-174, IProjectedCoordinateSystem branch).
    dist = np.linalg.norm(xy - centroid, axis=1)

    # MCPRange.cs:282: theDistances.Sort() - ascending order.
    order = np.argsort(dist, kind="stable")

    rows = []
    errors: list[str] = []
    for p in percentiles:
        try:
            # MCPRange.cs:253-254: need >=3 points at the desired percent level to draw a polygon.
            cutoff = int(np.floor(p * n))
            if cutoff < 3:
                raise ValueError(
                    "There are not enough points in the trajectory for the desired percent level - "
                    "need at least three points to draw an MCP polygon."
                )

            # MCPRange.cs:286-294: keep the `cutoff` closest-to-centroid fixes.
            kept = xy[order[:cutoff]]
            polygon = MultiPoint(kept).convex_hull  # MCPRange.cs:298-299: ITopologicalOperator2.ConvexHull()
            mcp_area = polygon.area  # MCPRange.cs:302-303: IArea.Area, raw native-CRS units (m^2 here)
            actual_percentile = 100.0 * cutoff / n  # MCPRange.cs:286

            rows.append(
                {
                    "ChosenPercentile": p * 100.0,  # MCPRange.cs:337
                    "ActualPercentile": actual_percentile,  # MCPRange.cs:336
                    "MCPArea": mcp_area,  # MCPRange.cs:338
                    "geometry": polygon,
                }
            )
        except Exception as ex:  # noqa: BLE001 - mirrors MCPRange.cs:344-348's catch-and-record-per-percentile
            msg = f"There was an error in the CalculateMCPRange function: {ex}"
            logger.warning(msg)
            errors.append(msg)

    gdf = gpd.GeoDataFrame(rows, crs=relocs_gdf.crs)
    return gdf, errors


def to_arcmet_schema(mcp_gdf: gpd.GeoDataFrame, subject_id: str, relocs_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Insert the identifier/date columns from MCPRange.cs's DefineFields() (lines 396-474),
    in the exact same field order: MovDataID, CalcID, StartDate, EndDate,
    ChosenPercentile, ActualPercentile, MCPArea.
    """
    out = mcp_gdf.copy()
    out.insert(0, "MovDataID", subject_id)
    out.insert(1, "CalcID", CALC_ID)
    # MCPRange.cs:320-335: non-moving-window analysis takes Start/End straight from the trajectory
    # (here: the raw relocations' own timestamp range - thePath.StartTime/EndTime are themselves
    # just the earliest/latest relocation timestamps).
    out.insert(2, "StartDate", relocs_gdf["fixtime"].min())
    out.insert(3, "EndDate", relocs_gdf["fixtime"].max())
    return out[["MovDataID", "CalcID", "StartDate", "EndDate", "ChosenPercentile", "ActualPercentile", "MCPArea", "geometry"]]


def build_summary_row(
    subject_id: str,
    relocs_gdf: gpd.GeoDataFrame,
    calc_initiated: dt.datetime,
    calc_finished: dt.datetime,
    result_fc_name: str,
    errors: list[str],
) -> dict:
    """Mirrors MCPRangeMethodResults.WriteCalcSummaryTable's fields (MCPRange.cs:34-98)."""
    return {
        "MovDataID": subject_id,
        "CalcID": CALC_ID,
        "CalcInitiated": calc_initiated,
        "CalcFinished": calc_finished,
        "StartDate": relocs_gdf["fixtime"].min(),
        "EndDate": relocs_gdf["fixtime"].max(),
        "ResultFCName": result_fc_name,
        "Errors": "None" if not errors else ";".join(errors),  # MCPRange.cs:51-64
    }


def run_subject(client, subject_name: str, subject_id: str, output_dir: str) -> dict:
    logger.info(f"=== {subject_name} (MCP) ===")
    calc_initiated = dt.datetime.now(dt.timezone.utc)  # MCPRange.cs:204: CalcInitiated = DateTime.UtcNow

    traj_gdf = bbmm.fetch_subject_trajectory(client, subject_id)
    relocs_gdf = bbmm.relocations_from_trajectory(traj_gdf)
    logger.info(f"{len(relocs_gdf)} relocations reconstructed from trajectory")

    mcp_gdf, errors = compute_mcp(relocs_gdf)

    slug = subject_name.lower().replace(" ", "_")
    os.makedirs(output_dir, exist_ok=True)
    result_fc_name = f"{slug}_MCP"

    arcmet_gdf = to_arcmet_schema(mcp_gdf, subject_id, relocs_gdf)
    csv_path = os.path.join(output_dir, f"{slug}_mcp_percentiles.csv")
    arcmet_gdf.drop(columns="geometry").to_csv(csv_path, index=False)
    logger.info(f"Wrote percentile table (ArcMET schema): {csv_path}")

    for _, row in mcp_gdf.sort_values("ChosenPercentile").iterrows():
        logger.info(
            f"  {row['ChosenPercentile']:>7.3f}% (actual {row['ActualPercentile']:.2f}%) "
            f"MCP area: {row['MCPArea'] / 1e6:.2f} km^2 ({row['MCPArea']:.0f} m^2)"
        )

    calc_finished = dt.datetime.now(dt.timezone.utc)  # MCPRange.cs:381: CalcFinished = DateTime.UtcNow
    summary_row = build_summary_row(subject_id, relocs_gdf, calc_initiated, calc_finished, result_fc_name, errors)
    summary_path = os.path.join(output_dir, f"{slug}_mcp_summary.csv")
    pd.DataFrame([summary_row]).to_csv(summary_path, index=False)
    logger.info(f"Wrote summary table: {summary_path}")

    # render_percentile_map expects "percentile"/"area_sqkm" columns (map-display
    # convenience only, not part of ArcMET's own schema/units); feed it the
    # pre-ArcMET-rename frame. rings=True converts these nested/cumulative hulls
    # into non-overlapping display rings, so opacity doesn't compound in the core.
    render_gdf = mcp_gdf.rename(columns={"ChosenPercentile": "percentile", "ActualPercentile": "actual_percentile"})
    render_gdf["area_sqkm"] = render_gdf["MCPArea"] / 1e6

    # Overlay the actual relocation points themselves (not a connecting line/Trajectory -
    # MCP never needed one; the raw fixes are what we actually have and computed the hull from).
    track_gdf = relocs_gdf[["geometry"]].copy()

    map_path = os.path.join(output_dir, f"{slug}_mcp_map.html")
    bbmm.render_percentile_map(
        render_gdf,
        map_path,
        colorize=True,
        filled=True,
        opacity=0.6,
        rings=True,
        track_gdf=track_gdf,
        track_color=[0, 0, 0],
        track_opacity=0.8,
        track_point_radius=1.5,
        track_on_top=True,
    )

    return summary_row


def main() -> None:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    for subject_name, subject_id in bbmm.SUBJECTS.items():
        run_subject(client, subject_name, subject_id, OUTPUT_DIR)


if __name__ == "__main__":
    main()
