# ETD Home Range

Computes a subject group's home range using the **Elliptical Time-Density (ETD)** model - a
trajectory-based, nonparametric estimate of the utilization distribution (UD), derived directly
from the group's movement behaviour (fixes → relocations → trajectory). This is the ETD slice of
the [subject-tracking workflow](https://github.com/ecoscope-platform-workflows-releases/subject-tracking),
stripped down to the home-range analysis plus a single percentile map - no dashboard.

## Outputs

Every run writes these artifacts to `$ECOSCOPE_WORKFLOWS_RESULTS`:

| File | Format | Contents |
|---|---|---|
| `etd_home_range_raster.tif` | GeoTIFF | The raw ETD utilization-distribution surface (continuous density, one band). |
| `etd_home_range_percentiles.parquet` / `.csv` | GeoParquet + CSV | One row per isopleth level: `percentile`, `geometry` (the percentile polygon), `area_sqkm`. |
| `etd_home_range_map.html` | HTML | Percentile choropleth map (RdYlGn by isopleth level) over base tile layers. |

The raster and the percentile table are computed from the same trajectory with the same grid
settings (CRS, cell size, nodata, max speed factor, expansion factor) - see the comment in
`spec.yaml` above the `Elliptical Time-Density (ETD)` task-group before changing one without the
other.

## Pipeline

```
Data Source → Time Range → Subject Group
  → Get Subject Observations (get_subjectgroup_observations)
  → Transform to Relocations (process_relocations)
  → Convert to Trajectory (relocations_to_trajectory)
  → Calculate ETD Percentiles (calculate_elliptical_time_density)   → GeoParquet + CSV
  → Generate ETD Raster (generate_etd_raster, ext-etd)              → GeoTIFF
  → Color + reproject to WGS84 + draw_map                           → HTML
```

`calculate_elliptical_time_density` is a built-in `ecoscope-platform` task and only returns the
percentile-area polygons - it doesn't expose the underlying raster. `generate_etd_raster` is a
small custom task (in [`wd-partner-tasks/src/ecoscope-workflows-ext-etd`](../wd-partner-tasks))
that calls the same `ecoscope.analysis.UD.calculate_etd_range` function with an `output_path`, so
the raw density surface can be persisted as GeoTIFF too.

The percentile geometry is computed in `EPSG:3857` (meters, for `area_sqkm`) but deck.gl expects
WGS84 lon/lat for its view state, so the map task-group reprojects a copy via `convert_crs` before
drawing - the analysis CRS and the display CRS are kept separate on purpose.

## Setup

`param.yaml` and `test-cases.yaml` are configured for the `mep` EarthRanger data source, subject
group `Salif`, over 2000-01-01 to 2020-01-01. To point this at a different data source, subject
group, or window, edit:

- `er_client.data_source.name` → your configured EarthRanger data source
- `subject_group_var.var` → the EarthRanger subject group to compute a home range for
- `time_range` → the analysis window

## Build & run

The `ecoscope-workflows-ext-etd` package isn't published yet, so `spec.yaml` points at a local
conda channel built from `wd-partner-tasks`.

```bash
make build     # builds ecoscope-workflows-ext-etd into /tmp/ecoscope-workflows-custom/release/artifacts
make compile   # compiles spec.yaml -> ecoscope-workflows-etd-workflow/
make run       # runs the workflow with param.yaml
```

Or step by step, without the local Makefile:

```bash
cd ../wd-partner-tasks
BUILD_PKG=etd pixi run build-release

cd ../ETD
pixi run compile-etd

cd ecoscope-workflows-etd-workflow
ECOSCOPE_WORKFLOWS_RESULTS="file:///tmp/workflows/etd/output" \
  pixi run ecoscope-workflows-etd-workflow run --config-file ../param.yaml --execution-mode sequential --no-mock-io
```

Once `ecoscope-workflows-ext-etd` is tagged (`git tag etd-v0.0.1 && git push --tags` in
`wd-partner-tasks`) and published to prefix.dev, switch the `requirements` entry in `spec.yaml`
from the local `file://` channel to `https://repo.prefix.dev/ecoscope-workflows-custom/` with a
real version, then recompile.
