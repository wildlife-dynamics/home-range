# ETD Home Range

Computes a subject group's home range using the **Elliptical Time-Density (ETD)** model - a
trajectory-based, nonparametric estimate of an animal's utilization distribution (UD), derived
directly from its own movement behaviour rather than a fitted statistical kernel. The method
builds "time-geography" ellipses between temporally adjacent GPS fixes, sized from a Weibull
distribution fit to the animal's own speed, and sums their overlap across the landscape to produce
a continuous time-density surface (Wall et al. 2014, *Methods in Ecology and Evolution*).

This workflow takes a subject group's fixes all the way from EarthRanger to a percentile-area
table, a percentile isopleth map, and the raw density-surface GeoTIFF - optionally split by
subject, time period, or spatial feature group - and wires the map and table into a dashboard.

## Outputs

Every run writes these to `$ECOSCOPE_WORKFLOWS_RESULTS`:

| Output | File | Contents |
|---|---|---|
| Dashboard | `result.json` | A map widget (percentile isopleth map) and a table widget (percentile/area table), with one view per group if groupers are used. |
| CSV | `etd_home_range_percentiles.csv` | One row per group per isopleth level: `percentile`, `geometry`, `area_sqkm`. |
| GeoTIFF | `etd_home_range_raster.tif` | The raw ETD utilization-distribution surface (continuous density, one band), one per group. |

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
