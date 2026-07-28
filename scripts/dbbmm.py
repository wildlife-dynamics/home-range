"""Dynamic Brownian Bridge Movement Model (dBBMM) home-range estimation.

Standalone prototype, extending bbmm.py's classic BBMM with the method of
Kranstauber, Kays, LaPoint, Wikelski & Safi (2012), "A dynamic Brownian
bridge movement model to estimate utilization distributions for
heterogeneous animal movement", J. Animal Ecology 81:738-746.

Why dynamic: classic BBMM assumes one Brownian motion variance (sigma_m^2)
for the *entire* track, which over- or under-estimates it in different
parts of a heterogeneous path (e.g. resting vs. travelling) and blurs the
resulting UD. dBBMM instead slides a window across the track and, at each
window, uses BIC to decide whether one sigma_m^2 fits the whole window
better than splitting it into two at some breakpoint (Kranstauber et al.
2012, Eqns 1-2). Estimates from all windows covering a given point are
averaged, giving a locally varying sigma_m^2 used in the same UD formula
as classic BBMM.

Reuses bbmm.py's data fetch, grid, and raster/percentile/map plumbing -
only the motion-variance estimation differs.
"""

from __future__ import annotations

import logging
import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import bbmm

logger = logging.getLogger(__name__)

# ── Configuration ────────────────────────────────────────────────────────────
# Window/margin must both be odd (Kranstauber et al. 2012): the odd-i-only
# likelihood needs an odd count of midpoints for a valid estimate, and margins
# guarantee >=3 non-breakpoint-eligible midpoints at each window edge.
WINDOW_SIZE = 31
MARGIN = 11
# The paper slides the window by one midpoint at a time. That's O(n) windows,
# each requiring O(window_size) breakpoint-candidate MLE fits - intractable in
# pure Python/scipy for tracks with tens of thousands of points. STRIDE trades
# some temporal resolution in the sigma_m^2(t) estimate for tractability; the
# window itself still covers every point, just with less overlap between
# neighbouring windows' averaged estimates.
STRIDE = 5

OUTPUT_DIR = os.environ.get("DBBMM_OUTPUT_DIR", "/tmp/workflows/dbbmm/output")


def _fit_sigma_m2(window: pd.DataFrame, location_error: float) -> tuple[float, float]:
    """MLE sigma_m^2 for a subset of midpoints; returns (sigma_m2, minimized -logL)."""
    result = minimize_scalar(
        bbmm._neg_log_likelihood,
        bounds=(1.0, 1_000_000.0),
        method="bounded",
        args=(window, location_error),
    )
    return float(result.x), float(result.fun)


def _bic(neg_log_lik: float, n_params: int, n_obs: int) -> float:
    return 2 * neg_log_lik + n_params * np.log(n_obs)


def estimate_dynamic_motion_variance(
    traj_gdf,
    location_error: float,
    window_size: int = WINDOW_SIZE,
    margin: int = MARGIN,
    stride: int = STRIDE,
) -> np.ndarray:
    """Per-segment sigma_m^2 array (Kranstauber et al. 2012's dBBMM), same length as
    the number of trajectory segments (one fewer than the number of points).
    """
    if window_size % 2 == 0 or margin % 2 == 0:
        raise ValueError("window_size and margin must both be odd (Kranstauber et al. 2012).")
    if margin * 2 >= window_size:
        raise ValueError("margin must be less than half of window_size.")

    midpoints = bbmm.compute_midpoints(traj_gdf)
    n = len(midpoints)
    n_segments = len(traj_gdf)  # each trajectory row is one segment; ground truth, not derived from midpoints

    if n < window_size:
        logger.warning(
            f"Track has only {n} valid midpoints (< window_size={window_size}); "
            "falling back to a single classic-BBMM estimate for the whole track."
        )
        sigma_m2, _ = _fit_sigma_m2(midpoints, location_error)
        return np.full(n_segments, sigma_m2)

    sums = np.zeros(n)
    counts = np.zeros(n)

    starts = list(range(0, n - window_size + 1, stride))
    if starts[-1] != n - window_size:
        starts.append(n - window_size)  # make sure the tail of the track is covered too

    candidate_offsets = [k for k in range(margin, window_size - margin) if k % 2 == 1]

    logger.info(f"dBBMM: sliding {len(starts)} windows (size={window_size}, margin={margin}, stride={stride}) over {n} midpoints")
    for w_i, w_start in enumerate(starts):
        window = midpoints.iloc[w_start : w_start + window_size]
        sigma_whole, neg_ll_whole = _fit_sigma_m2(window, location_error)
        bic_whole = _bic(neg_ll_whole, n_params=1, n_obs=window_size)

        best_bic_split = np.inf
        best_split = None
        for offset in candidate_offsets:
            left, right = window.iloc[:offset], window.iloc[offset:]
            if len(left) < 3 or len(right) < 3:
                continue
            sigma_left, neg_ll_left = _fit_sigma_m2(left, location_error)
            sigma_right, neg_ll_right = _fit_sigma_m2(right, location_error)
            bic_split = _bic(neg_ll_left + neg_ll_right, n_params=2, n_obs=window_size)
            if bic_split < best_bic_split:
                best_bic_split = bic_split
                best_split = (offset, sigma_left, sigma_right)

        core = range(margin, window_size - margin)
        if best_split is not None and best_bic_split < bic_whole:
            offset, sigma_left, sigma_right = best_split
            for k in core:
                idx = w_start + k
                sums[idx] += sigma_left if k < offset else sigma_right
                counts[idx] += 1
        else:
            for k in core:
                idx = w_start + k
                sums[idx] += sigma_whole
                counts[idx] += 1

        if (w_i + 1) % 200 == 0:
            logger.info(f"  ...processed {w_i + 1}/{len(starts)} windows")

    per_midpoint_sigma = pd.Series(np.where(counts > 0, sums / np.maximum(counts, 1), np.nan))
    n_uncovered = per_midpoint_sigma.isna().sum()
    if n_uncovered:
        logger.info(f"  {n_uncovered} midpoint(s) at track edges uncovered by any window core; forward/backward-filling")
    per_midpoint_sigma = per_midpoint_sigma.ffill().bfill()

    segment_sigma = pd.Series(np.full(n_segments, np.nan))
    segment_sigma.loc[midpoints["segment_before"].to_numpy()] = per_midpoint_sigma.to_numpy()
    segment_sigma.loc[midpoints["segment_after"].to_numpy()] = per_midpoint_sigma.to_numpy()
    segment_sigma = segment_sigma.ffill().bfill()  # segments between two odd-i midpoints (even-i gaps)

    logger.info(
        f"dBBMM sigma_m^2: min={segment_sigma.min():.2f}, median={segment_sigma.median():.2f}, "
        f"max={segment_sigma.max():.2f} m^2/s (vs one classic-BBMM value for the whole track)"
    )
    return segment_sigma.to_numpy()


def run_subject(client, subject_name: str, subject_id: str, output_dir: str) -> None:
    logger.info(f"=== {subject_name} (dBBMM) ===")
    traj_gdf = bbmm.fetch_subject_trajectory(client, subject_id)
    logger.info(f"{len(traj_gdf)} trajectory segments fetched")

    sigma_m2_per_segment = estimate_dynamic_motion_variance(traj_gdf, bbmm.LOCATION_ERROR_METERS)

    from ecoscope.analysis.UD import grid_size_from_geographic_extent

    pixel_size = grid_size_from_geographic_extent(traj_gdf, scale_factor=bbmm.GRID_SCALE_FACTOR)
    logger.info(f"Grid pixel size: {pixel_size} m")

    raster_data, raster_profile = bbmm.compute_bbmm_ud(
        traj_gdf, sigma_m2_per_segment, bbmm.LOCATION_ERROR_METERS, pixel_size
    )

    slug = subject_name.lower().replace(" ", "_")
    os.makedirs(output_dir, exist_ok=True)
    raster_path = os.path.join(output_dir, f"{slug}_dbbmm_raster.tif")

    write_array = raster_data.data.copy()
    write_array[np.isnan(write_array) | (write_array == 0)] = raster_profile.nodata_value
    bbmm.raster.RasterPy.write(write_array, fp=raster_path, **raster_profile)
    logger.info(f"Wrote GeoTIFF: {raster_path}")

    percentiles_gdf = bbmm.get_percentile_area(bbmm.PERCENTILES, raster_data, subject_id=subject_name)
    percentiles_gdf = percentiles_gdf.set_geometry("geometry", crs=raster_profile.crs)
    percentiles_gdf["area_sqkm"] = percentiles_gdf.geometry.area / 1e6

    csv_path = os.path.join(output_dir, f"{slug}_dbbmm_percentiles.csv")
    percentiles_gdf.drop(columns="geometry").to_csv(csv_path, index=False)
    logger.info(f"Wrote percentile table: {csv_path}")

    for _, row in percentiles_gdf.sort_values("percentile").iterrows():
        logger.info(f"  {row['percentile']:>7.3f}% isopleth area: {row['area_sqkm']:.2f} km^2")

    map_path = os.path.join(output_dir, f"{slug}_dbbmm_map.html")
    bbmm.render_percentile_map(percentiles_gdf, map_path)


def main() -> None:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    for subject_name, subject_id in bbmm.SUBJECTS.items():
        run_subject(client, subject_name, subject_id, OUTPUT_DIR)


if __name__ == "__main__":
    main()
