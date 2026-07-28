"""Ad-hoc comparison: rings-correction on vs off, for both ETD and MCP, on the same
subject (Habiba) - to decide whether rings-correction should be a shared, unified
toggle across both methods in the ETD/MCP home-range selector, or left MCP-only.
"""

from __future__ import annotations

import logging
import os

import bbmm
import home_range_dryrun as hrd
from ecoscope.trajectory import Trajectory

logger = logging.getLogger(__name__)

OUTPUT_DIR = os.environ.get("RINGS_COMPARISON_OUTPUT_DIR", "/tmp/workflows/rings_comparison/output")


def main() -> list[str]:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    subject_name, subject_id = "Habiba", bbmm.SUBJECTS["Habiba"]

    traj_gdf = bbmm.fetch_subject_trajectory(client, subject_id)
    relocations_gdf = Trajectory(gdf=traj_gdf).to_relocations().gdf

    etd_result = hrd.run_etd(traj_gdf)
    mcp_result = hrd.run_mcp(relocations_gdf)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    paths = []

    # (opacity_label, etd_opacity, mcp_opacity) - 0.25 is the original opacity where the
    # nested-layer compounding problem was first discovered and looked dramatic; 0.6/0.7
    # are what we've actually settled on for MCP/ETD, where the effect is much subtler.
    opacity_variants = [("low_opacity", 0.25, 0.25), ("chosen_opacity", 0.7, 0.6)]

    for opacity_label, etd_opacity, mcp_opacity in opacity_variants:
        for rings in (True, False):
            suffix = "rings_on" if rings else "rings_off"

            etd_path = os.path.join(OUTPUT_DIR, f"habiba_etd_{opacity_label}_{suffix}.html")
            bbmm.render_percentile_map(
                etd_result, etd_path, colorize=True, filled=True, stroked=False, opacity=etd_opacity, rings=rings
            )
            logger.info(f"Wrote: {etd_path}")
            paths.append(etd_path)

            mcp_path = os.path.join(OUTPUT_DIR, f"habiba_mcp_{opacity_label}_{suffix}.html")
            bbmm.render_percentile_map(
                mcp_result,
                mcp_path,
                colorize=True,
                filled=True,
                opacity=mcp_opacity,
                rings=rings,
                track_gdf=relocations_gdf[["geometry"]],
                track_color=[0, 0, 0],
                track_opacity=0.8,
                track_point_radius=1.5,
                track_on_top=True,
                legend_title="Home Range %",
            )
            logger.info(f"Wrote: {mcp_path}")
            paths.append(mcp_path)

    return paths


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    for path in main():
        print(path)
