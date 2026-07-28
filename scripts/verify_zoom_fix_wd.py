"""Verifies the _ASSUMED_VIEWPORT_WIDTH_PX change (1200 -> 650) in wd's
compute_fitted_view_state (ecoscope_workflows_ext_wd/tasks/_home_range.py)
against Habiba's real trajectory - the same production task used by
spec.yaml's "Compute Fitted Map View per Group" step.

Prints the resulting zoom for both the old (1200px) and new (650px) viewport
assumptions side by side, and renders both maps to HTML so the difference can
be seen visually, not just as a number.
"""

from __future__ import annotations

import logging

import bbmm
import math
import pydeck as pdk

from ecoscope.platform.tasks.results import create_geojson_layer, draw_map
from ecoscope.platform.tasks.results._map_utils import TileLayer
from ecoscope.platform.tasks.results._pydeck import GeoJSONLayerStyle, LegendFromDataframe, LegendStyle
from ecoscope.platform.tasks.transformation import apply_color_map, convert_column_values_to_string, convert_crs, map_columns

from ecoscope_workflows_ext_wd.tasks._home_range import EtdMethodArgs, call_home_range_from_args, compute_fitted_view_state

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

SUBJECT_NAME = "Habiba"
SUBJECT_ID = bbmm.SUBJECTS[SUBJECT_NAME]
_TILE_WIDTH_PX = 256


def render(view_state: pdk.ViewState, title: str, output_path: str, df) -> None:
    layer = create_geojson_layer(
        geodataframe=df,
        layer_style=GeoJSONLayerStyle(filled=True, stroked=False, get_fill_color="percentile_color", opacity=0.7),
        legend=LegendFromDataframe(
            title="Home Range %", label_column="Percentile", label_suffix=" %",
            color_column="percentile_color", sort="ascending",
        ),
        tooltip_columns=["Percentile", "Area"],
        zoom=False,
    )
    html = draw_map(
        geo_layers=[layer],
        tile_layers=[TileLayer(layer_name="TERRAIN", opacity=1.0)],
        static=False,
        output_type="html",
        view_state=view_state,
        max_zoom=15,
        title=title,
        legend_style=LegendStyle(placement="bottom-right"),
    )
    with open(output_path, "w") as f:
        f.write(html)
    logger.info(f"wrote {output_path}  ->  {view_state}")


def main() -> None:
    client = bbmm.EarthRangerConnection.client_from_named_connection(bbmm.DATA_SOURCE)
    traj_gdf = bbmm.fetch_subject_trajectory(client, SUBJECT_ID)
    logger.info(f"{SUBJECT_NAME}: {len(traj_gdf)} trajectory segments fetched")

    args = EtdMethodArgs(
        etd_opacity=0.7,
        etd_crs="EPSG:3857",
        max_speed_factor=1.05,
        etd_expansion_factor=1.3,
        etd_percentiles=[50.0, 60.0, 70.0, 80.0, 90.0, 99.999],
    )
    result = call_home_range_from_args(traj_gdf, args)
    logger.info(f"{SUBJECT_NAME}: computed {len(result)} ETD percentile polygons")

    df = convert_column_values_to_string(result, columns=["percentile"])
    df = apply_color_map(df, input_column_name="percentile", colormap="RdYlGn_r", output_column_name="percentile_color")
    df = convert_crs(df, crs="EPSG:4326")
    df = map_columns(
        df, drop_columns=[], retain_columns=[],
        rename_columns={"percentile": "Percentile", "area_sqkm": "Area"}, raise_if_not_found=True,
    )

    layer_def = create_geojson_layer(
        geodataframe=df,
        layer_style=GeoJSONLayerStyle(filled=True, stroked=False, get_fill_color="percentile_color", opacity=0.7),
        legend=None,
        tooltip_columns=None,
        zoom=False,
        data_url=None,
    )

    # AFTER: production compute_fitted_view_state with the new 650px assumption.
    vs_650 = compute_fitted_view_state(geo_layers=layer_def, max_zoom=15)

    # OLD: same bbox math, but with the previous 1200px assumption.
    bounds = df.total_bounds
    bbox = [[bounds[0], bounds[1]], [bounds[2], bounds[3]]]
    naive_zoom = pdk.data_utils.viewport_helpers.bbox_to_zoom_level(bbox)
    zoom_1200 = min(15, naive_zoom + math.log2(1200 / _TILE_WIDTH_PX))
    vs_1200 = pdk.ViewState(longitude=vs_650.longitude, latitude=vs_650.latitude, zoom=zoom_1200)

    print(f"naive (256px, no correction):     zoom={naive_zoom:.2f}")
    print(f"OLD assumption (1200px viewport): zoom={zoom_1200:.2f}")
    print(f"NEW assumption (650px viewport):  zoom={vs_650.zoom:.2f}")
    print(f"difference: {zoom_1200 - vs_650.zoom:.2f} fewer zoom levels with the new assumption")

    render(pdk.ViewState(longitude=vs_650.longitude, latitude=vs_650.latitude, zoom=zoom_1200),
           "OLD (1200px assumption - over-zoomed for a ~650px widget)",
           "/tmp/habiba_zoomfix_1200.html", df)
    render(vs_650, "NEW (650px assumption)", "/tmp/habiba_zoomfix_650.html", df)


if __name__ == "__main__":
    main()
