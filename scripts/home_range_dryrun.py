"""Dry run of the actual ETD/MCP home-range selector now implemented in ecoscope
(ecoscope.platform.tasks.config: EtdMethodArgs, McpMethodArgs, set_home_range_args,
call_home_range_from_combined_params, plus the rings/legend/stroked/points-overlay
getters) - exercising the real selector pathway end to end against real EarthRanger
data for Habiba and Salif Keita, before wiring it into ETD/spec.yaml itself.

Reuses bbmm.py's data fetch and map plumbing; the home-range computations and all
method-selection logic come straight from ecoscope's platform task layer.
"""

from __future__ import annotations

import logging
import os

import bbmm
from ecoscope.platform.tasks.config import (
    EtdMethodArgs,
    McpMethodArgs,
    call_home_range_from_combined_params,
    get_legend_title_from_combined_params,
    get_opacity_from_combined_params,
    get_rings_correction_from_combined_params,
    get_stroked_from_combined_params,
    relocations_for_points_overlay,
    set_home_range_args,
)
from ecoscope.platform.tasks.preprocessing import convert_trajectory_to_relocations

logger = logging.getLogger(__name__)

PERCENTILES = bbmm.PERCENTILES  # full set: 50/60/70/80/90/95/99.999

OUTPUT_DIR = os.environ.get("HOME_RANGE_DRYRUN_OUTPUT_DIR", "/tmp/workflows/home_range_dryrun/output")


def log_results(subject_name: str, method: str, result) -> None:
    for _, row in result.sort_values("percentile").iterrows():
        logger.info(f"  [{subject_name}] {method} {row['percentile']:>6.2f}%  area: {row['area_sqkm']:.2f} km^2")


def render(result, points_overlay, combined_params, subject_name: str, method_label: str, output_dir: str) -> str:
    """Every rendering choice (stroked, opacity, legend title, whether points show at
    all) is read from the selector's own combined_params, exactly as the spec.yaml
    wiring will - nothing here is hardcoded per method by this script."""
    slug = subject_name.lower().replace(" ", "_")
    map_path = os.path.join(output_dir, f"{slug}_{method_label}_dryrun_map.html")
    bbmm.render_percentile_map(
        result,
        map_path,
        colorize=True,
        filled=True,
        stroked=get_stroked_from_combined_params(combined_params),
        opacity=get_opacity_from_combined_params(combined_params),
        rings=get_rings_correction_from_combined_params(combined_params),
        track_gdf=points_overlay if len(points_overlay) else None,
        track_color=[0, 0, 0],
        track_opacity=0.8,
        track_point_radius=1.5,
        track_on_top=True,
        legend_title=get_legend_title_from_combined_params(combined_params),
    )
    logger.info(f"Wrote map: {map_path}")
    return map_path


def run_subject(client, subject_name: str, subject_id: str, output_dir: str) -> list[str]:
    logger.info(f"=== {subject_name} ===")
    traj_gdf = bbmm.fetch_subject_trajectory(client, subject_id)
    logger.info(f"{len(traj_gdf)} trajectory segments fetched")

    relocations_gdf = convert_trajectory_to_relocations(traj_gdf)
    logger.info(f"{len(relocations_gdf)} relocations derived via convert_trajectory_to_relocations")

    etd_combined = set_home_range_args(opacity=0.7, method_args=EtdMethodArgs(percentiles=PERCENTILES))
    mcp_combined = set_home_range_args(opacity=0.6, method_args=McpMethodArgs(crs=bbmm.CRS, percentiles=PERCENTILES))

    etd_result = call_home_range_from_combined_params(traj_gdf, relocations_gdf, etd_combined)
    mcp_result = call_home_range_from_combined_params(traj_gdf, relocations_gdf, mcp_combined)
    log_results(subject_name, "ETD", etd_result)
    log_results(subject_name, "MCP", mcp_result)

    etd_columns, mcp_columns = set(etd_result.columns), set(mcp_result.columns)
    logger.info(f"[{subject_name}] ETD columns: {sorted(etd_columns)}, MCP columns: {sorted(mcp_columns)}")
    assert etd_columns == mcp_columns, "shared downstream chain would break!"

    etd_points = relocations_for_points_overlay(relocations_gdf, etd_combined)
    mcp_points = relocations_for_points_overlay(relocations_gdf, mcp_combined)
    logger.info(f"[{subject_name}] points overlay rows - ETD: {len(etd_points)}, MCP: {len(mcp_points)}")
    assert len(etd_points) == 0
    assert len(mcp_points) == len(relocations_gdf)

    os.makedirs(output_dir, exist_ok=True)
    etd_map = render(etd_result, etd_points, etd_combined, subject_name, "etd", output_dir)
    mcp_map = render(mcp_result, mcp_points, mcp_combined, subject_name, "mcp", output_dir)
    return [etd_map, mcp_map]


def main() -> list[str]:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    paths = []
    for subject_name, subject_id in bbmm.SUBJECTS.items():
        paths.extend(run_subject(client, subject_name, subject_id, OUTPUT_DIR))
    return paths


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    for path in main():
        print(path)
