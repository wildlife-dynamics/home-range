# Home Range

Computes a subject group's **home range** (the area an animal traverses during its normal
activities of foraging, mating, and caring for young, per Burt, 1943) from its EarthRanger tracking
data.

Choose one of three estimation methods per run:

- **Elliptical Time-Density (ETD)** — a trajectory-based, nonparametric estimate of an animal's
  utilization distribution (UD), derived directly from its own movement behaviour rather than a
  fitted statistical kernel. Builds "time-geography" ellipses between temporally adjacent GPS
  fixes, sized from a Weibull distribution fit to the animal's own speed, and sums their overlap
  across the landscape to produce a continuous time-density surface (Wall et al. 2014, *Methods in
  Ecology and Evolution*).
- **Minimum Convex Polygon (MCP)** — the classic, purely geometric estimator (Mohr 1947): ranks
  fixes by distance from the centroid and draws the convex hull enclosing the closest fraction at
  each percentile level, dropping distant excursions as outliers.
- **Brownian Bridge Movement Model (BBMM)** — a UD estimate like ETD, but modeled as a Brownian
  bridge between consecutive fixes rather than a fixed-size ellipse, so its uncertainty scales with
  the time gap between fixes (Horne et al. 2007, *Ecology*).

Only the selected method's own settings are shown in the form; ETD and BBMM also produce a GeoTIFF
of the underlying utilization-distribution surface, since MCP has no equivalent (it's pure
geometry, not a density surface).

## Outputs

Every run generates a dashboard and a percentiles dataframe; ETD and BBMM additionally produce a
GeoTIFF of the home range's utilization-distribution surface.

## Setup

Edit `param.yaml`:

- `er_client.data_source.name` → your EarthRanger data source
- `subject_observations.subject_group_name` → the subject group to compute a home range for
- `time_range` → the analysis window
- `home_range_args.args` → the estimation method (ETD/MCP/BBMM) and its settings

## Build & run

```bash
pixi run compile-etd

cd ecoscope-workflows-home-range-workflow
ECOSCOPE_WORKFLOWS_RESULTS="file:///tmp/workflows/etd/output" \
  pixi run ecoscope-workflows-home-range-workflow run --config-file ../param.yaml --execution-mode sequential --no-mock-io
```
