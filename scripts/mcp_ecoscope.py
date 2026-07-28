"""Run ecoscope's own calculate_mcp_range (ecoscope.analysis.UD, feature/mcp-home-range
branch of /Users/zak/Documents/w-dynamics/ecoscope) against real EarthRanger data for
Habiba and Salif Keita - a demo/smoke test of the library implementation, as opposed to
mcp.py's standalone ArcMET-faithful port.

Reuses bbmm.py's data fetch and map plumbing; only the MCP computation itself comes from
ecoscope now.
"""

from __future__ import annotations

import logging
import os

import bbmm
from ecoscope.analysis.UD import calculate_mcp_range

logger = logging.getLogger(__name__)

PERCENTILES = [50.0, 80.0, 99.0]

OUTPUT_DIR = os.environ.get("MCP_ECOSCOPE_OUTPUT_DIR", "/tmp/workflows/mcp_ecoscope/output")


def run_subject(client, subject_name: str, subject_id: str, output_dir: str) -> str:
    logger.info(f"=== {subject_name} (ecoscope calculate_mcp_range) ===")

    traj_gdf = bbmm.fetch_subject_trajectory(client, subject_id)
    relocs_gdf = bbmm.relocations_from_trajectory(traj_gdf)
    logger.info(f"{len(relocs_gdf)} relocations reconstructed from trajectory")

    result = calculate_mcp_range(
        relocations=relocs_gdf,
        percentile_levels=PERCENTILES,
        crs=bbmm.CRS,
        subject_id=subject_name,
    )
    result["area_sqkm"] = result.area / 1_000_000.0

    slug = subject_name.lower().replace(" ", "_")
    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, f"{slug}_ecoscope_mcp.csv")
    result.drop(columns="geometry").to_csv(csv_path, index=False)
    logger.info(f"Wrote: {csv_path}")

    for _, row in result.iterrows():
        logger.info(
            f"  {row['percentile']:>6.2f}% (actual {row['actual_percentile']:.2f}%) area: {row['area_sqkm']:.2f} km^2"
        )

    map_path = os.path.join(output_dir, f"{slug}_ecoscope_mcp_map.html")
    bbmm.render_percentile_map(
        result,
        map_path,
        colorize=True,
        filled=True,
        opacity=0.6,
        rings=True,
        track_gdf=relocs_gdf[["geometry"]],
        track_color=[0, 0, 0],
        track_opacity=0.8,
        track_point_radius=1.5,
        track_on_top=True,
        legend_title="Home Range %",
    )
    logger.info(f"Wrote map: {map_path}")
    return map_path


def main() -> list[str]:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    return [run_subject(client, subject_name, subject_id, OUTPUT_DIR) for subject_name, subject_id in bbmm.SUBJECTS.items()]


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    for path in main():
        print(path)
