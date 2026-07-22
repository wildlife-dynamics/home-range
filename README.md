# Home Range

Computes a subject group's **home range** (the area an animal traverses during its normal
activities of foraging, mating, and caring for young, per Burt, 1943) from its EarthRanger tracking
data.

Currently implemented using the **Elliptical Time-Density (ETD)** model: a trajectory-based,
nonparametric estimate of an animal's utilization distribution (UD), derived directly from its own
movement behaviour rather than a fitted statistical kernel. The method builds "time-geography"
ellipses between temporally adjacent GPS fixes, sized from a Weibull distribution fit to the
animal's own speed, and sums their overlap across the landscape to produce a continuous
time-density surface (Wall et al. 2014, *Methods in Ecology and Evolution*).

We're hoping to support additional home-range estimation methods (e.g. kernel density estimation,
minimum convex polygon) alongside ETD going forward.

## Outputs

Every run generates a dashboard, a percentiles dataframe, and a GeoTIFF of the home range.

## Setup

Edit `param.yaml`:

- `er_client.data_source.name` → your EarthRanger data source
- `subject_observations.subject_group_name` → the subject group to compute a home range for
- `time_range` → the analysis window

## Build & run

```bash
pixi run compile-etd

cd ecoscope-workflows-etd-workflow
ECOSCOPE_WORKFLOWS_RESULTS="file:///tmp/workflows/etd/output" \
  pixi run ecoscope-workflows-etd-workflow run --config-file ../param.yaml --execution-mode sequential --no-mock-io
```
