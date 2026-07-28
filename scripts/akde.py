"""Autocorrelated Kernel Density Estimation (aKDE) home-range estimation.

Standalone prototype, implementing the stationary-process AKDE derived in
Fleming, Fagan, Mueller, Olson, Leimgruber & Calabrese (2015), "Rigorous
home range estimation with movement data: a new autocorrelated kernel
density estimator", Ecology 96(5):1182-1188 (Appendix B.2.1 - the
stationary-process special case of their general derivation).

Why autocorrelated: conventional KDE treats every fix as an independent
sample, which badly underestimates home range size for realistically
autocorrelated GPS data - nearby-in-time fixes are nearly duplicates, not
independent draws, so denser sampling makes a conventional KDE's home
range *smaller* (and more wrong), not more accurate. AKDE instead:

  1. estimates the population's stationary covariance sigma_0 and its
     empirical semi-variance function (SVF) gamma(tau) directly from the
     trajectory (Eqn B.33/B.36);
  2. fits a smooth parametric autocorrelation function to that empirical
     SVF - here an Ornstein-Uhlenbeck exponential decay, g(tau) = 1 -
     exp(-tau/tau_p), the simplest and most common choice;
  3. numerically minimizes the mean integrated square error (MISE, Eqn
     B.40) over a single scalar bandwidth multiplier h, correcting the
     kernel bandwidth upward for whatever autocorrelation the fitted ACF
     implies;
  4. produces a standard (but now correctly-sized) weighted KDE using the
     corrected bandwidth sigma_B = h^2 * sigma_0.

This is Fleming et al.'s own foundational, stationary-process special
case, not the full `ctmm` R package's more general moving/non-stationary
OUF-process + Kalman-filter machinery - a faithful but simpler rendition
of the same core idea, well-suited to a single home-range estimate per
individual.

Reuses bbmm.py's data fetch, grid, and raster/percentile/map plumbing -
only the bandwidth estimation and final density formula differ (no
bridging between fixes at all; every point contributes one fixed-width
kernel).
"""

from __future__ import annotations

import logging
import os

import numpy as np
from scipy.optimize import minimize_scalar

import bbmm
from ecoscope.io import raster

logger = logging.getLogger(__name__)

# ── Configuration ────────────────────────────────────────────────────────────
MAX_VARIOGRAM_PAIRS = 200_000  # random-subsample point pairs, for tractability on large tracks (Salif Keita: ~2.4e8 total pairs)
MAX_LAG_QUANTILE = 0.5  # fit the ACF only to lags at/below this quantile of all sampled pairwise lags (short lags are what pin down tau_p)
N_LAG_BINS = 30

OUTPUT_DIR = os.environ.get("AKDE_OUTPUT_DIR", "/tmp/workflows/akde/output")


def _sample_pairs(n_points: int, max_pairs: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Distinct-index (i, j) pairs, exhaustive if small enough, otherwise a random subsample."""
    total_pairs = n_points * (n_points - 1) // 2
    if total_pairs <= max_pairs:
        return np.triu_indices(n_points, k=1)
    i_idx = rng.integers(0, n_points, size=max_pairs)
    j_idx = rng.integers(0, n_points, size=max_pairs)
    valid = i_idx != j_idx
    i_idx, j_idx = i_idx[valid], j_idx[valid]
    return np.minimum(i_idx, j_idx), np.maximum(i_idx, j_idx)


def estimate_akde_bandwidth(traj_gdf, rng: np.random.Generator | None = None) -> dict:
    """Fit the stationary AKDE bandwidth (Fleming et al. 2015, Appendix B.2.1).

    Returns sigma_0 (2x2 stationary covariance), tau_p (OU position
    autocorrelation timescale, seconds), h (optimal scalar bandwidth
    multiplier, Eqn B.40), and sigma_B = h^2 * sigma_0 (corrected kernel
    covariance used for the final KDE).
    """
    rng = rng or np.random.default_rng(0)
    xy, t = bbmm._extract_points(traj_gdf)
    n = len(xy)

    mu = xy.mean(axis=0)
    centered = xy - mu
    sigma_0 = (centered.T @ centered) / n
    trace_sigma_0 = np.trace(sigma_0)

    i_idx, j_idx = _sample_pairs(n, MAX_VARIOGRAM_PAIRS, rng)
    lag = np.abs(t[i_idx] - t[j_idx])
    sq_dist = np.sum((centered[i_idx] - centered[j_idx]) ** 2, axis=1)
    # trace(2*gamma(tau)) = E[|r(t)-r(t')|^2] (Eqn B.33/B.36); with gamma(tau) = g(tau)*sigma_0
    # (Eqn B.39, constant-anisotropy assumption), g(tau) = mean(sq_dist) / (2 * trace(sigma_0)).
    g_empirical = sq_dist / (2 * trace_sigma_0)

    max_lag = np.quantile(lag, MAX_LAG_QUANTILE)
    keep = lag <= max_lag
    lag_fit, g_fit = lag[keep], g_empirical[keep]

    bin_edges = np.linspace(0, max_lag, N_LAG_BINS + 1)
    bin_idx = np.clip(np.digitize(lag_fit, bin_edges) - 1, 0, N_LAG_BINS - 1)
    binned_g = np.array(
        [g_fit[bin_idx == b].mean() if (bin_idx == b).any() else np.nan for b in range(N_LAG_BINS)]
    )
    binned_tau = (bin_edges[:-1] + bin_edges[1:]) / 2
    valid_bins = ~np.isnan(binned_g)
    binned_tau, binned_g = binned_tau[valid_bins], binned_g[valid_bins]

    def _sse(log_tau_p: float) -> float:
        tau_p = np.exp(log_tau_p)
        pred = 1 - np.exp(-binned_tau / tau_p)
        return float(np.sum((binned_g - pred) ** 2))

    tau_result = minimize_scalar(_sse, bounds=(np.log(1.0), np.log(max(t[-1] - t[0], 2.0))), method="bounded")
    tau_p = float(np.exp(tau_result.x))
    logger.info(f"Fitted OU position autocorrelation timescale (tau_p) = {tau_p:.1f} s ({tau_p / 3600:.2f} h)")

    q = 2  # spatial dimension
    g_at_lag = 1 - np.exp(-lag / tau_p)

    def _mise(h: float) -> float:
        h2 = h * h
        term1 = np.mean(1.0 / (2 * g_at_lag + 2 * h2) ** (q / 2))
        term2 = 2.0 / (2 + h2) ** (q / 2)
        term3 = 1.0 / 2 ** (q / 2)
        return float(term1 - term2 + term3)

    h_result = minimize_scalar(_mise, bounds=(1e-3, 10.0), method="bounded")
    h = float(h_result.x)
    logger.info(f"AKDE bandwidth multiplier h = {h:.4f} (sigma_B = h^2 * sigma_0)")

    sigma_B = (h**2) * sigma_0
    return {
        "mu": mu,
        "sigma_0": sigma_0,
        "tau_p": tau_p,
        "h": h,
        "sigma_B": sigma_B,
        "n_pairs_used": len(lag),
    }


def compute_akde_ud(
    traj_gdf,
    sigma_B: np.ndarray,
    pixel_size: float,
    expansion_factor: float = bbmm.EXPANSION_FACTOR,
):
    """Weighted KDE: sum one fixed-covariance Gaussian kernel per fix, no bridging between fixes."""
    raster_profile, col_centers, row_centers = bbmm._build_grid(traj_gdf, pixel_size, expansion_factor)
    num_rows, num_cols = raster_profile.rows, raster_profile.columns
    ud = np.zeros((num_rows, num_cols), dtype=np.float64)

    xy, _ = bbmm._extract_points(traj_gdf)
    n = len(xy)

    sigma_inv = np.linalg.inv(sigma_B)
    det_sigma_B = np.linalg.det(sigma_B)
    norm_const = 1.0 / (2 * np.pi * np.sqrt(det_sigma_B))

    max_sigma = np.sqrt(np.linalg.eigvalsh(sigma_B).max())
    pad = bbmm.MAX_BRIDGE_RADIUS_SIGMA * max_sigma + pixel_size

    for i in range(n):
        x0, y0 = xy[i]
        col_mask = (col_centers >= x0 - pad) & (col_centers <= x0 + pad)
        row_mask = (row_centers >= y0 - pad) & (row_centers <= y0 + pad)
        if not col_mask.any() or not row_mask.any():
            continue

        local_x, local_y = col_centers[col_mask], row_centers[row_mask]
        X, Y = np.meshgrid(local_x, local_y)
        dx, dy = X - x0, Y - y0
        quad = dx * dx * sigma_inv[0, 0] + 2 * dx * dy * sigma_inv[0, 1] + dy * dy * sigma_inv[1, 1]
        density = norm_const * np.exp(-0.5 * quad)

        row_idx, col_idx = np.where(row_mask)[0], np.where(col_mask)[0]
        ud[np.ix_(row_idx, col_idx)] += density

        if (i + 1) % 5000 == 0:
            logger.info(f"  ...processed {i + 1}/{n} points")

    ud /= n  # (1/n) sum over point-kernels

    cell_area = pixel_size * pixel_size
    total_mass = ud.sum() * cell_area
    if total_mass > 0:
        ud /= total_mass

    return raster.RasterData(data=ud.astype("float32"), crs=raster_profile.crs, transform=raster_profile.transform), raster_profile


def run_subject(client, subject_name: str, subject_id: str, output_dir: str) -> None:
    logger.info(f"=== {subject_name} (aKDE) ===")
    traj_gdf = bbmm.fetch_subject_trajectory(client, subject_id)
    logger.info(f"{len(traj_gdf)} trajectory segments fetched")

    fit = estimate_akde_bandwidth(traj_gdf)
    logger.info(
        f"sigma_0 (stationary covariance, m^2): {fit['sigma_0'].tolist()} | "
        f"h={fit['h']:.4f} | tau_p={fit['tau_p']:.1f}s | pairs used={fit['n_pairs_used']}"
    )

    from ecoscope.analysis.UD import grid_size_from_geographic_extent

    pixel_size = grid_size_from_geographic_extent(traj_gdf, scale_factor=bbmm.GRID_SCALE_FACTOR)
    logger.info(f"Grid pixel size: {pixel_size} m")

    raster_data, raster_profile = compute_akde_ud(traj_gdf, fit["sigma_B"], pixel_size)

    slug = subject_name.lower().replace(" ", "_")
    os.makedirs(output_dir, exist_ok=True)
    raster_path = os.path.join(output_dir, f"{slug}_akde_raster.tif")

    write_array = raster_data.data.copy()
    write_array[np.isnan(write_array) | (write_array == 0)] = raster_profile.nodata_value
    raster.RasterPy.write(write_array, fp=raster_path, **raster_profile)
    logger.info(f"Wrote GeoTIFF: {raster_path}")

    percentiles_gdf = bbmm.get_percentile_area(bbmm.PERCENTILES, raster_data, subject_id=subject_name)
    percentiles_gdf = percentiles_gdf.set_geometry("geometry", crs=raster_profile.crs)
    percentiles_gdf["area_sqkm"] = percentiles_gdf.geometry.area / 1e6

    csv_path = os.path.join(output_dir, f"{slug}_akde_percentiles.csv")
    percentiles_gdf.drop(columns="geometry").to_csv(csv_path, index=False)
    logger.info(f"Wrote percentile table: {csv_path}")

    for _, row in percentiles_gdf.sort_values("percentile").iterrows():
        logger.info(f"  {row['percentile']:>7.3f}% isopleth area: {row['area_sqkm']:.2f} km^2")

    map_path = os.path.join(output_dir, f"{slug}_akde_map.html")
    bbmm.render_percentile_map(percentiles_gdf, map_path)


def main() -> None:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    for subject_name, subject_id in bbmm.SUBJECTS.items():
        run_subject(client, subject_name, subject_id, OUTPUT_DIR)


if __name__ == "__main__":
    main()
