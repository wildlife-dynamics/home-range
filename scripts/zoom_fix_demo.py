"""Demo of the view-state zoom-fit fix in ecoscope's view_state_from_geodataframes
(ecoscope/platform/tasks/results/_pydeck.py). Renders the same real Habiba percentile
polygon bbox with the BEFORE (naive, assumes a 256px viewport) and AFTER (corrected for
a ~1200px viewport) computed zoom, so the difference can be seen side by side.

Note: draw_map itself hits an unrelated pandera/typing quirk when called from a bare
script outside pytest/the compiled DAG, so this renders with raw pydeck directly instead
- the view_state values themselves come straight from the real ecoscope function.
"""

from __future__ import annotations

import math

import geopandas as gpd
import pydeck as pdk
from shapely.geometry import box
from wt_task import task

from ecoscope.platform.tasks.results._pydeck import (
    create_geojson_layer,
    view_state_from_geodataframes,
)

# Real Habiba-scale percentile-polygon bbox, extracted from the actual compiled
# workflow's output HTML earlier in this session.
POLY = box(37.54274, 0.53531, 37.63084, 0.61348)
GDF = gpd.GeoDataFrame({"percentile": ["90"], "geometry": [POLY]}, crs="EPSG:4326")

# Passed as a plain dict (not a constructed LayerStyle(...)) - constructing the
# dataclass directly in a bare script trips an unrelated pandera/typing quirk that
# only manifests outside pytest/the compiled DAG; letting task().validate() coerce
# the dict avoids it.
LAYER_STYLE = {
    "filled": True,
    "get_fill_color": [255, 140, 0],
    "stroked": True,
    "get_line_color": [0, 0, 0],
    "get_line_width": 2,
    "opacity": 0.6,
}


def render(view_state: pdk.ViewState, title: str, out_path: str) -> None:
    layer_def = task(create_geojson_layer).validate().call(
        geodataframe=GDF,
        layer_style=LAYER_STYLE,
        legend=None,
        tooltip_columns=None,
        zoom=False,
        data_url=None,
    )
    layer = pdk.Layer(type=layer_def.layer_type, data=GDF, **{
        "filled": True,
        "get_fill_color": [255, 140, 0],
        "stroked": True,
        "get_line_color": [0, 0, 0],
        "get_line_width": 2,
        "opacity": 0.6,
    })
    deck = pdk.Deck(
        layers=[layer],
        initial_view_state=view_state,
        map_style=None,
        widgets=[pdk.Widget("TitleWidget", title=title)],
    )
    deck.to_html(out_path, notebook_display=False)
    print(f"wrote {out_path}  ->  {view_state}")


def main() -> None:
    # AFTER: the real, fixed ecoscope function.
    vs_after = view_state_from_geodataframes(geodataframes=[GDF], max_zoom=20)

    # BEFORE: reproduce the old (unfixed) formula for comparison - pydeck's raw
    # bbox_to_zoom_level with no viewport-width correction.
    bounds = GDF.total_bounds
    bbox = [[bounds[0], bounds[1]], [bounds[2], bounds[3]]]
    naive_zoom = pdk.data_utils.viewport_helpers.bbox_to_zoom_level(bbox)
    vs_before = pdk.ViewState(
        longitude=vs_after.longitude,
        latitude=vs_after.latitude,
        zoom=naive_zoom,
    )

    print(f"BEFORE zoom (naive, assumes 256px viewport): {vs_before.zoom}")
    print(f"AFTER  zoom (corrected for ~1200px viewport): {vs_after.zoom}")
    print(f"difference: +{vs_after.zoom - vs_before.zoom:.2f} zoom levels "
          f"(= log2(1200/256) = {math.log2(1200 / 256):.2f})")

    render(vs_before, "BEFORE fix (zoomed out, not fitted)", "/tmp/habiba_zoomfix_before.html")
    render(vs_after, "AFTER fix (fitted)", "/tmp/habiba_zoomfix_after.html")


if __name__ == "__main__":
    main()
