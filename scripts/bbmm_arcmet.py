"""Brownian Bridge Movement Model (BBMM), ArcMET-comparable port.

A line-for-line port of ArcMET's own C# implementation at
/Users/zak/Documents/w-dynamics/arcmet/ArcMET_Base/ArcMET_Base/BBMMRange.cs
(class BBMMRange.CalculateBBMMRange + BBLikelihood + GoldenSectionSearch),
extending bbmm.py (which implements the *textbook* Horne et al. 2007 method,
validated against R's CRAN `BBMM` package and Kranstauber et al. 2012 - no
ArcMET involved at all).

Unlike mcp.py (a pure bug-for-bug ArcMET replica, because ArcMET's own
MCPRange.cs turned out to have no formula bugs to begin with - "textbook"
and "ArcMET" are simply the same algorithm for MCP), BBMMRange.cs contains
two genuine bugs, plus several behavioral choices (gap handling, probability
trimming, normalization convention, discretization granularity) that are
implementation details rather than the core Horne et al. formula. So every
point of divergence below is its own parameter, defaulting to the *textbook*
(bbmm.py-equivalent) value - call with `**ARCMET_FIDELITY_PRESET` to instead
reproduce ArcMET's own behavior exactly, bugs included.

Two genuine bugs in BBMMRange.cs (confirmed by direct comparison, not just
its own comments):

  1. BBLikelihood's own bridge-variance formula (BBMMRange.cs:611):
         bbVariance = T*alpha*(1-alpha)*mobilityVariance
                      + (alpha*alpha + (1 - alpha*alpha)) * locVariance
     The parenthesized term algebraically simplifies to exactly 1 for *any*
     alpha (alpha^2 + (1-alpha^2) = 1) - a copy-paste typo that silently
     collapses the location-error weighting to a constant, instead of the
     correct (1-alpha)^2 + alpha^2 used two lines away in the raster loop
     itself (BBMMRange.cs:435). Only affects the MLE fit of sigma_m^2, not
     the final raster.
     Toggle: `reproduce_arcmet_likelihood_bug` (default False).

  2. The main raster loop's density term (BBMMRange.cs:444):
         theta = (1 / Math.Sqrt(2*pi*bbVariance)) * exp(-distanceSquared/(2*bbVariance))
     `1/sqrt(2*pi*var)` is the 1-D Gaussian normalization constant, applied
     to a 2-D (x,y) distance - the correct 2-D bivariate-normal constant is
     `1/(2*pi*var)`, which is what BBLikelihood itself correctly uses two
     functions away (BBMMRange.cs:612). This one distorts the actual UD
     raster ArcMET writes to disk, not just the MLE fit.
     Toggle: `reproduce_arcmet_density_normalization_bug` (default False).

Behavioral/config divergences (not bugs - deliberate ArcMET choices, each
its own parameter, textbook default first / ArcMET value second):

  - `max_data_gap_seconds` (None / 14400.0, BBMMRange.cs:21,424): ArcMET
    excludes segments with a gap >= this from BOTH the density accumulation
    and the totalMinuteSum normalization denominator (BBMMRange.cs:424-451).
    Textbook default keeps every segment regardless of gap size.
  - `probability_cutoff` (0.0 / 1e-6, BBMMRange.cs:23,483-500): ArcMET zeroes
    (NoData) any final cell below this, then renormalizes the survivors to
    sum to 1. Textbook default keeps the full continuous tail.
  - `weight_by_time_step` (True / False, BBMMRange.cs:449): ArcMET's inner
    loop never multiplies theta by the time step before accumulating - it
    sums raw instantaneous densities, then forces the whole grid to
    renormalize to 1 as a compensating hack (BBMMRange.cs:455-474), which
    implicitly overweights segments with more sub-steps. Textbook default
    Riemann-integrates each segment's own time contribution (`* dt`) before
    accumulating, so segment duration doesn't distort relative weighting.
  - `normalize_by_cell_area` (True / False, BBMMRange.cs:456-474): ArcMET's
    "divide by totalMinuteSum, rescale to sum to 1" normalizes the RAW cell
    values to sum to 1, ignoring cell area entirely, so the raster's absolute
    values are not resolution-independent as a continuous density. Textbook
    default additionally divides by cell area, giving a proper continuous PDF.
    NOTE: this only changes the raw GeoTIFF's absolute values - it has NO
    effect on the percentile-area table, since get_percentile_area
    (ecoscope/analysis/percentile.py) works entirely off cumulative-sum
    *ratios* of sorted pixel values, which are invariant to a single global
    rescaling like dividing every cell by the same cell_area constant.
  - `time_step_seconds` (60.0 / 600.0, BBMMRange.cs:20 `mTimeStep = 10`
    minutes): discretization granularity for the bridge-integration loop.
  - `max_sd` (6.0 / 5.0, BBMMRange.cs:22 `mMaxSD = 5`) together with
    `apply_circular_truncation` (False / True, BBMMRange.cs:436,442): ArcMET
    zeroes any cell farther than `max_sd` standard deviations from the exact
    per-timestep bridge mean (`maxDist = MaxSD^2 * bbVariance`), recomputed
    every time-step. Textbook default has no such per-cell check at all -
    bbmm.py's own compute_bbmm_ud only ever restricts to a padded rectangular
    window per segment (sized off the segment's own max variance across all
    its sub-steps) and accumulates every cell in it unconditionally,
    including the window's corners (which sit farther out than the circle
    max_sd describes). With max_sd generous enough this is numerically
    negligible either way (Gaussian tails beyond 5-6 sigma are ~1e-8 of
    peak density), but `apply_circular_truncation=False` is kept as the
    default specifically so the textbook path stays bit-for-bit identical
    to bbmm.py's own output, not just "close".
  - `use_arcmet_golden_section_search` (False / True, BBMMRange.cs:637-816):
    ArcMET fits sigma_m^2 via a hand-rolled translation of Numerical
    Recipes' bracket + golden-section-search routines (Press et al. 3rd
    ed., pg. 491, 495), with ax=1/bx=1,000,000 (in per-minute units -
    rescaled to per-second units here to match this port's second-based
    time convention throughout) and tol=3^-8 (BBMMRange.cs:642 - the
    in-code comment claims this should be sqrt(machine epsilon), but
    3^-8 =~ 1.5e-4 is not that; kept verbatim as ArcMET's own value, not
    "corrected", since replicating it is the whole point of this toggle).
    Textbook default uses scipy.optimize.minimize_scalar (Brent's method),
    a more robust, well-tested bounded optimizer - both should converge to
    the same minimum of whatever objective function they're given, so this
    only matters if you want the optimizer *process* itself to match, not
    just the input/output formulas.

Time units: this port keeps everything in SECONDS throughout (matching
bbmm.py, and _extract_points's own seconds-based timestamps) rather than
ArcMET's minutes - this is inconsequential to results (sigma_m^2 is a free
MLE-fit parameter; using seconds vs minutes just rescales its fitted value
by a constant factor, as long as the same convention is used consistently
in both the MLE-fit step and the raster step, which it is here) so it is
NOT exposed as a toggle. Only `use_arcmet_golden_section_search`'s internal
bracket bounds are rescaled (by /60) to remain equivalent to ArcMET's own
minute-based bounds despite the second-based convention used here.

ArcMET's own raster-extent/grid setup (RasterDataHelper.SetupOutputFileGeodatabaseRaster,
not present in BBMMRange.cs itself, defined elsewhere in ArcMET_Base and not
available to inspect) can't be replicated byte-for-byte - this reuses
bbmm.py's own `_build_grid` (bbox + expansion_factor + pixel_size), the same
convention already used for ETD/BBMM/MCP elsewhere in this repo. Likewise,
BBMMRange.cs contains no percentile-area table computation at all (that
appears to live in a shared moving-window/raster percentile helper not in
this file) - percentile-area output here reuses ecoscope's own
get_percentile_area, exactly as bbmm.py's textbook version does.

Reuses bbmm.py's data fetch, grid, and raster/percentile/map plumbing -
only the motion-variance estimation and UD computation differ.
"""

from __future__ import annotations

import logging
import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import bbmm
from ecoscope.analysis.percentile import get_percentile_area
from ecoscope.io import raster

logger = logging.getLogger(__name__)

# ── Configuration ────────────────────────────────────────────────────────────

OUTPUT_DIR = os.environ.get("BBMM_ARCMET_OUTPUT_DIR", "/tmp/workflows/bbmm_arcmet/output")

# Every one of ArcMET's values in one place, for `run_subject(..., **ARCMET_FIDELITY_PRESET)`.
ARCMET_FIDELITY_PRESET = dict(
    reproduce_arcmet_likelihood_bug=True,
    reproduce_arcmet_density_normalization_bug=True,
    max_data_gap_seconds=14400.0,  # BBMMRange.cs:21 default (4 hours)
    probability_cutoff=1e-6,  # BBMMRange.cs:23 default
    weight_by_time_step=False,
    normalize_by_cell_area=False,
    time_step_seconds=600.0,  # BBMMRange.cs:20 default (10 minutes)
    # ArcMET checks maxDist against the FULL raster grid every time-step (no
    # windowing concept at all in its own code) - window_padding_sigma must be
    # >= max_sd here so this port's vectorized local window doesn't clip cells
    # the exact circular check (below) would otherwise have accepted.
    window_padding_sigma=5.0,
    max_sd=5.0,  # BBMMRange.cs:22 default
    apply_circular_truncation=True,
    use_arcmet_golden_section_search=True,
)


# ── Brownian motion variance (MLE) ──────────────────────────────────────────


def _arcmet_bridge_variance_for_likelihood(jump_time: np.ndarray, alpha: np.ndarray, sigma_m2: float, location_variance: float) -> np.ndarray:
    """BBMMRange.cs:611, verbatim - including its own bug. `(alpha**2 + (1 - alpha**2))`
    algebraically equals 1 for any alpha; written out in full (not simplified to
    `+ location_variance`) so the bug stays visually traceable to that exact line."""
    return jump_time * alpha * (1 - alpha) * sigma_m2 + (alpha * alpha + (1 - alpha * alpha)) * location_variance


def _arcmet_neg_log_likelihood(sigma_m2: float, midpoints: pd.DataFrame, location_variance: float, reproduce_arcmet_likelihood_bug: bool) -> float:
    """BBMMRange.cs:589-626 (BBLikelihood), verbatim including its exact per-point
    try/except-style density==0 skip and its correct 2-D density formula (only
    the *variance* formula is buggy, not the density formula built from it)."""
    jump_time = midpoints["jump_time"].to_numpy()
    alpha = midpoints["alpha"].to_numpy()
    d2 = midpoints["d2"].to_numpy()

    if reproduce_arcmet_likelihood_bug:
        variance = _arcmet_bridge_variance_for_likelihood(jump_time, alpha, sigma_m2, location_variance)
    else:
        variance = bbmm._bridge_variance(jump_time, alpha, sigma_m2, np.sqrt(location_variance))

    # BBMMRange.cs:612: density = (1/(2*pi*bbVariance)) * exp(-0.5*(distanceSquared/bbVariance)) - correct 2-D form.
    density = (1.0 / (2 * np.pi * variance)) * np.exp(-0.5 * d2 / variance)
    # BBMMRange.cs:613-617: only accumulate where density != 0 (guards log(0)).
    nonzero = density != 0
    log_density = np.log(density[nonzero])
    sum_density = log_density.sum()
    return -1.0 * sum_density  # BBMMRange.cs:619: likelihood = -1 * SumDensity


class _GoldenSectionSearch:
    """Literal translation of BBMMRange.cs's GoldenSectionSearch class
    (BBMMRange.cs:637-816) - Numerical Recipes' Bracket + golden-section
    Minimize routines (Press et al. 3rd ed., pg. 491/495)."""

    TOL = 3.0**-8  # BBMMRange.cs:642, kept verbatim (not "corrected" to sqrt(machine epsilon))
    GOLD = 1.618034
    GLIMIT = 100.0
    TINY = 1.0e-20
    R = 0.61803399
    C = 1 - R

    def __init__(self, objective):
        self._f = objective  # objective(x) -> float, all other args already bound

    def bracket(self, ax: float, bx: float) -> tuple[float, float, float]:
        """BBMMRange.cs:700-781."""
        f = self._f
        fa, fb = f(ax), f(bx)
        if fb > fa:
            ax, bx = bx, ax
            fa, fb = fb, fa
        cx = bx + self.GOLD * (bx - ax)
        if cx < 0:
            cx = bx - bx * 0.99  # BBMMRange.cs:727, "From Horne et al. 2007"
        fc = f(cx)
        while fb > fc:
            r = (bx - ax) * (fb - fc)
            q = (bx - cx) * (fb - fa)
            denom = 2 * self._sign(max(abs(q - r), self.TINY), q - r)
            u = bx - ((bx - cx) * q - (bx - ax) * r) / denom
            ulim = bx + self.GLIMIT * (cx - bx)
            if (bx - u) * (u - cx) > 0:
                fu = f(u)
                if fu < fc:
                    return bx, u, cx
                elif fu > fb:
                    return ax, bx, u
                u = cx + self.GOLD * (cx - bx)
                fu = f(u)
            elif (cx - u) * (u - ulim) > 0:
                fu = f(u)
                if fu < fc:
                    bx, cx, u = cx, u, u + self.GOLD * (u - cx)
                    fb, fc, fu = fc, fu, f(u)
            elif (u - ulim) * (ulim - cx) >= 0:
                u = ulim
                fu = f(u)
            else:
                u = cx + self.GOLD * (cx - bx)
                fu = f(u)
            ax, bx, cx = bx, cx, u
            fa, fb, fc = fb, fc, fu
        return ax, bx, cx

    def minimize(self, ax: float, bx: float, cx: float) -> float:
        """BBMMRange.cs:644-698."""
        f = self._f
        x0, x3 = ax, cx
        if abs(cx - bx) > abs(bx - ax):
            x1, x2 = bx, bx + self.C * (cx - bx)
        else:
            x2, x1 = bx, bx - self.C * (bx - ax)
        f1, f2 = f(x1), f(x2)
        while abs(x3 - x0) > self.TOL * (abs(x1) + abs(x2)):
            if f2 < f1:
                x0, x1, x2 = x1, x2, self.R * x2 + self.C * x3
                f1, f2 = f2, f(x2)
            else:
                x3, x2, x1 = x2, x1, self.R * x1 + self.C * x0
                f2, f1 = f1, f(x1)
        return x1 if f1 < f2 else x2

    @staticmethod
    def _sign(a: float, b: float) -> float:
        return abs(a) if b >= 0 else -abs(a)


def estimate_motion_variance(
    traj_gdf,
    location_error: float,
    reproduce_arcmet_likelihood_bug: bool = False,
    use_arcmet_golden_section_search: bool = False,
) -> float:
    """MLE for a single, whole-track sigma_m^2. Reuses bbmm.py's own odd-index
    leave-one-out midpoint selection (compute_midpoints) - that scheme matches
    BBMMRange.cs's `for (i=1; i<numOfLocs-1; i+=2)` exactly regardless of which
    variance/optimizer variant is used."""
    midpoints = bbmm.compute_midpoints(traj_gdf)
    if midpoints.empty:
        raise ValueError("Not enough interior fixes with valid time gaps to estimate motion variance.")
    location_variance = location_error**2

    if not use_arcmet_golden_section_search:
        result = minimize_scalar(
            _arcmet_neg_log_likelihood,
            bounds=(1.0, 1_000_000.0),
            method="bounded",
            args=(midpoints, location_variance, reproduce_arcmet_likelihood_bug),
        )
        sigma_m2 = float(result.x)
    else:
        # BBMMRange.cs:309-313: ax=1, bx=1,000,000 in per-MINUTE units (Horne et al.'s
        # own R code); this port is second-based throughout, so /60 keeps the same
        # physical bracket (bbVariance ~ time * sigma_m2, so an equivalent sigma_m2
        # bound scales by 1/60 when the time unit is 60x smaller).
        ax, bx = 1.0 / 60.0, 1_000_000.0 / 60.0

        def objective(x: float) -> float:
            return _arcmet_neg_log_likelihood(x, midpoints, location_variance, reproduce_arcmet_likelihood_bug)

        search = _GoldenSectionSearch(objective)
        ax_b, bx_b, cx_b = search.bracket(ax, bx)
        sigma_m2 = search.minimize(ax_b, bx_b, cx_b)

    logger.info(
        f"Estimated Brownian motion variance (sigma_m^2) = {sigma_m2:.2f} m^2/s "
        f"(from {len(midpoints)} fixes, arcmet_likelihood_bug={reproduce_arcmet_likelihood_bug}, "
        f"arcmet_search={use_arcmet_golden_section_search})"
    )
    return sigma_m2


# ── UD raster computation ───────────────────────────────────────────────────


def compute_bbmm_ud(
    traj_gdf,
    sigma_m2: float,
    location_error: float,
    pixel_size: float,
    time_step_seconds: float = 60.0,
    expansion_factor: float = bbmm.EXPANSION_FACTOR,
    max_data_gap_seconds: float | None = None,
    probability_cutoff: float = 0.0,
    weight_by_time_step: bool = True,
    normalize_by_cell_area: bool = True,
    reproduce_arcmet_density_normalization_bug: bool = False,
    window_padding_sigma: float = bbmm.MAX_BRIDGE_RADIUS_SIGMA,
    max_sd: float = 5.0,
    apply_circular_truncation: bool = False,
    max_steps_per_segment: int = bbmm.MAX_STEPS_PER_SEGMENT,
):
    """BBMMRange.cs:259-521 (CalculateBBMMRange)'s main raster loop, parameterized
    per-divergence - see module docstring. All defaults are the textbook
    (bbmm.py-equivalent) values; pass **ARCMET_FIDELITY_PRESET for ArcMET's own."""
    raster_profile, col_centers, row_centers = bbmm._build_grid(traj_gdf, pixel_size, expansion_factor)
    num_rows, num_cols = raster_profile.rows, raster_profile.columns
    ud = np.zeros((num_rows, num_cols), dtype=np.float64)

    xy, t = bbmm._extract_points(traj_gdf)
    location_variance = location_error**2

    n_segments = len(xy) - 1
    total_time_sum = 0.0  # BBMMRange.cs:400 totalMinuteSum (named generically here: seconds, not minutes)

    for i in range(n_segments):
        z0, z1 = xy[i], xy[i + 1]
        time_lag = t[i + 1] - t[i]
        if time_lag <= 0:
            continue
        # BBMMRange.cs:424: segments with a gap >= MaxDataGapSeconds are excluded
        # entirely - from both density accumulation AND the normalization sum.
        if max_data_gap_seconds is not None and time_lag >= max_data_gap_seconds:
            continue

        total_time_sum += time_lag

        n_steps = min(max_steps_per_segment, max(2, int(np.ceil(time_lag / time_step_seconds)) + 1))
        alphas = np.linspace(0.0, 1.0, n_steps)
        dt = time_lag / (n_steps - 1)

        variances = bbmm._bridge_variance(time_lag, alphas, sigma_m2, location_error)
        max_sigma = np.sqrt(variances.max())
        # window_padding_sigma sizes the padded rectangular window (bbmm.py's own
        # MAX_BRIDGE_RADIUS_SIGMA concept - a vectorization/performance device, not
        # part of ArcMET's algorithm at all). max_sd (ArcMET's own MaxSD) is a
        # SEPARATE, only-if-apply_circular_truncation per-cell radius check below -
        # conflating the two was the earlier bug in this file (window defaulted to
        # max_sd's value instead of bbmm.py's own 4.0, silently widening the window
        # and picking up extra far-tail cells even in the "textbook" default path).
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
            dist_sq = (X - mu_x) ** 2 + (Y - mu_y) ** 2
            if reproduce_arcmet_density_normalization_bug:
                # BBMMRange.cs:444, verbatim bug: 1-D normalization constant on a 2-D distance.
                theta = (1.0 / np.sqrt(2 * np.pi * variance)) * np.exp(-dist_sq / (2 * variance))
            else:
                theta = (1.0 / (2 * np.pi * variance)) * np.exp(-dist_sq / (2 * variance))
            if apply_circular_truncation:
                # BBMMRange.cs:436,442: exact per-timestep circular truncation (maxDist =
                # MaxSD^2 * bbVariance). Textbook default (bbmm.py) has no such per-cell
                # check at all - every cell in the segment's padded rectangular window
                # gets a contribution, including the window's corners (further out than
                # this circle). Opt-in only, so the textbook path stays bit-for-bit
                # identical to bbmm.py's own compute_bbmm_ud.
                theta = np.where(dist_sq <= (max_sd**2 * variance), theta, 0.0)
            local_accum += theta

        if weight_by_time_step:
            local_accum = local_accum * dt  # textbook: Riemann-sum the time integral over this segment
        # else: BBMMRange.cs:449 - accumulated as raw per-tm-step density, uncorrected for step width.

        row_idx = np.where(row_mask)[0]
        col_idx = np.where(col_mask)[0]
        ud[np.ix_(row_idx, col_idx)] += local_accum

        if (i + 1) % 5000 == 0:
            logger.info(f"  ...processed {i + 1}/{n_segments} segments")

    if total_time_sum <= 0:
        raise ValueError("No segments survived max_data_gap_seconds filtering - nothing to compute a UD from.")

    # BBMMRange.cs:455-465: divide by totalMinuteSum (here: total_time_sum), then
    # rescale so the grid sums to exactly 1 - ArcMET's own two-step normalization,
    # applied here regardless of weight_by_time_step/normalize_by_cell_area so the
    # probability_cutoff step below always operates on a properly-scaled [0,1] grid.
    ud = ud / total_time_sum
    grid_sum = ud.sum()
    if grid_sum > 0:
        ud = ud / grid_sum

    if probability_cutoff > 0:
        # BBMMRange.cs:477-500: zero out (NoData) any cell below the cutoff, then
        # renormalize the survivors to sum to 1.
        survivors = ud >= probability_cutoff
        new_sum = ud[survivors].sum()
        ud = np.where(survivors, ud, 0.0)
        if new_sum > 0:
            ud = ud / new_sum

    if normalize_by_cell_area:
        # textbook: convert the [sums-to-1 discrete PMF] into a proper continuous
        # PDF (integrates to 1 over the raster's real area), by dividing by cell
        # area - ArcMET's own normalization (above) never does this.
        cell_area = pixel_size * pixel_size
        ud = ud / cell_area

    return raster.RasterData(data=ud.astype("float32"), crs=raster_profile.crs, transform=raster_profile.transform), raster_profile


# ── Main ─────────────────────────────────────────────────────────────────────


def run_subject(
    client,
    subject_name: str,
    subject_id: str,
    output_dir: str,
    location_error_meters: float = bbmm.LOCATION_ERROR_METERS,
    **arcmet_kwargs,
) -> None:
    """`**arcmet_kwargs` forwards to both estimate_motion_variance and
    compute_bbmm_ud - pass nothing for the textbook (bbmm.py-equivalent)
    defaults, or **ARCMET_FIDELITY_PRESET for ArcMET's own values."""
    variance_kwargs = {
        k: arcmet_kwargs[k]
        for k in ("reproduce_arcmet_likelihood_bug", "use_arcmet_golden_section_search")
        if k in arcmet_kwargs
    }
    ud_kwargs = {k: v for k, v in arcmet_kwargs.items() if k not in variance_kwargs}

    logger.info(f"=== {subject_name} (BBMM, ArcMET-comparable port) ===")
    traj_gdf = bbmm.fetch_subject_trajectory(client, subject_id)
    logger.info(f"{len(traj_gdf)} trajectory segments fetched")

    sigma_m2 = estimate_motion_variance(traj_gdf, location_error_meters, **variance_kwargs)

    from ecoscope.analysis.UD import grid_size_from_geographic_extent

    pixel_size = grid_size_from_geographic_extent(traj_gdf, scale_factor=bbmm.GRID_SCALE_FACTOR)
    logger.info(f"Grid pixel size: {pixel_size} m")

    raster_data, raster_profile = compute_bbmm_ud(traj_gdf, sigma_m2, location_error_meters, pixel_size, **ud_kwargs)

    slug = subject_name.lower().replace(" ", "_")
    os.makedirs(output_dir, exist_ok=True)
    raster_path = os.path.join(output_dir, f"{slug}_bbmm_arcmet_raster.tif")

    write_array = raster_data.data.copy()
    write_array[np.isnan(write_array) | (write_array == 0)] = raster_profile.nodata_value
    raster.RasterPy.write(write_array, fp=raster_path, **raster_profile)
    logger.info(f"Wrote GeoTIFF: {raster_path}")

    percentiles_gdf = get_percentile_area(bbmm.PERCENTILES, raster_data, subject_id=subject_name)
    percentiles_gdf = percentiles_gdf.set_geometry("geometry", crs=raster_profile.crs)
    percentiles_gdf["area_sqkm"] = percentiles_gdf.geometry.area / 1e6

    csv_path = os.path.join(output_dir, f"{slug}_bbmm_arcmet_percentiles.csv")
    percentiles_gdf.drop(columns="geometry").to_csv(csv_path, index=False)
    logger.info(f"Wrote percentile table: {csv_path}")

    for _, row in percentiles_gdf.sort_values("percentile").iterrows():
        logger.info(f"  {row['percentile']:>7.3f}% isopleth area: {row['area_sqkm']:.2f} km^2")

    map_path = os.path.join(output_dir, f"{slug}_bbmm_arcmet_map.html")
    # stroked=False: BBMM is a density-surface method like ETD, not a hull method
    # like MCP - many small isopleth fragments look like a messy quilted mesh with
    # an outline on every one (matches ETD's own production etd_pct_layer config).
    bbmm.render_percentile_map(percentiles_gdf, map_path, stroked=False)


def main() -> None:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    for subject_name, subject_id in bbmm.SUBJECTS.items():
        run_subject(client, subject_name, subject_id, OUTPUT_DIR)


if __name__ == "__main__":
    main()
