# Run with the local spatialdata installation from the repo root:
#   asv run --config benchmarks/asv.conf.json --environment existing:python3.13 HEAD^!
#
# Or compare specific commits:
#   asv run --config benchmarks/asv.conf.json HEAD^!
#
# Show results:
#   asv show --config benchmarks/asv.conf.json HEAD
#
# HTML report:
#   asv publish --config benchmarks/asv.conf.json
#   asv preview --config benchmarks/asv.conf.json

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, box

import spatialdata as sd
from spatialdata import bounding_box_query, polygon_query
from spatialdata.models import PointsModel, ShapesModel


def _make_points(n: int, rng: np.random.Generator) -> sd.SpatialData:
    xy = rng.uniform(0, 1000, (n, 2))
    return PointsModel.parse(pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}))


def _make_circles(n: int, rng: np.random.Generator) -> gpd.GeoDataFrame:
    xy = rng.uniform(0, 1000, (n, 2))
    radii = rng.uniform(5, 20, n)
    return ShapesModel.parse(
        gpd.GeoDataFrame(
            {"geometry": [Point(x, y) for x, y in xy], "radius": radii}
        )
    )


def _make_query_boxes(n_boxes: int) -> tuple[np.ndarray, np.ndarray]:
    """Evenly tile the [0, 1000]^2 space into n_boxes non-overlapping boxes."""
    step = 1000.0 / n_boxes
    mins = np.array([[i * step, i * step] for i in range(n_boxes)])
    maxs = mins + step * 0.8
    return mins, maxs


class SingleBoxQueryPoints:
    """
    Condition 1 — single rectangular query on **points**.

    Compares:
      * ``bounding_box_query`` (native implementation)
      * ``polygon_query`` with a box-shaped polygon (same spatial extent)
    """

    params = [1_000, 50_000, 500_000]
    param_names = ["n_points"]

    timeout = 120

    def setup(self, n_points):
        rng = np.random.default_rng(0)
        self.points = _make_points(n_points, rng)
        self.min_coord = [200.0, 200.0]
        self.max_coord = [700.0, 700.0]
        self.box_polygon = box(200, 200, 700, 700)

    def time_bounding_box_query(self, n_points):
        bounding_box_query(
            self.points,
            axes=("x", "y"),
            min_coordinate=self.min_coord,
            max_coordinate=self.max_coord,
            target_coordinate_system="global",
        ).compute()

    def time_polygon_query(self, n_points):
        polygon_query(
            self.points,
            polygon=self.box_polygon,
            target_coordinate_system="global",
        ).compute()


class SingleBoxQueryShapes:
    """
    Condition 1 — single rectangular query on **circles (shapes)**.

    Compares:
      * ``bounding_box_query`` (native implementation)
      * ``polygon_query`` with a box-shaped polygon (same spatial extent)
    """

    params = [500, 5_000, 50_000]
    param_names = ["n_shapes"]

    timeout = 120

    def setup(self, n_shapes):
        rng = np.random.default_rng(0)
        self.circles = _make_circles(n_shapes, rng)
        self.min_coord = [200.0, 200.0]
        self.max_coord = [700.0, 700.0]
        self.box_polygon = box(200, 200, 700, 700)

    def time_bounding_box_query(self, n_shapes):
        bounding_box_query(
            self.circles,
            axes=("x", "y"),
            min_coordinate=self.min_coord,
            max_coordinate=self.max_coord,
            target_coordinate_system="global",
        )

    def time_polygon_query(self, n_shapes):
        polygon_query(
            self.circles,
            polygon=self.box_polygon,
            target_coordinate_system="global",
        )


class MultiBboxQueryPoints:
    """
    Condition 2 — multiple rectangular queries on **points**.

    Compares:
      * ``bounding_box_query`` with a (n_boxes, 2) coordinate array (vectorised)
      * ``polygon_query`` called once per box in a loop
    """

    params = ([50_000], [5, 20, 50])
    param_names = ["n_points", "n_boxes"]

    timeout = 300

    def setup(self, n_points, n_boxes):
        rng = np.random.default_rng(0)
        self.points = _make_points(n_points, rng)
        self.mins, self.maxs = _make_query_boxes(n_boxes)
        self.box_polygons = [
            box(self.mins[i][0], self.mins[i][1], self.maxs[i][0], self.maxs[i][1])
            for i in range(n_boxes)
        ]

    def time_bounding_box_query_vectorized(self, n_points, n_boxes):
        results = bounding_box_query(
            self.points,
            axes=("x", "y"),
            min_coordinate=self.mins,
            max_coordinate=self.maxs,
            target_coordinate_system="global",
        )
        for r in results:
            r.compute()

    def time_polygon_query_loop(self, n_points, n_boxes):
        for poly in self.box_polygons:
            polygon_query(
                self.points,
                polygon=poly,
                target_coordinate_system="global",
            ).compute()


class MultiBboxQueryShapes:
    """
    Condition 2 — multiple rectangular queries on **circles (shapes)**.

    Compares:
      * ``bounding_box_query`` with a (n_boxes, 2) coordinate array (vectorised)
      * ``polygon_query`` called once per box in a loop
    """

    params = ([5_000], [5, 20, 50])
    param_names = ["n_shapes", "n_boxes"]

    timeout = 300

    def setup(self, n_shapes, n_boxes):
        rng = np.random.default_rng(0)
        self.circles = _make_circles(n_shapes, rng)
        self.mins, self.maxs = _make_query_boxes(n_boxes)
        self.box_polygons = [
            box(self.mins[i][0], self.mins[i][1], self.maxs[i][0], self.maxs[i][1])
            for i in range(n_boxes)
        ]

    def time_bounding_box_query_vectorized(self, n_shapes, n_boxes):
        bounding_box_query(
            self.circles,
            axes=("x", "y"),
            min_coordinate=self.mins,
            max_coordinate=self.maxs,
            target_coordinate_system="global",
        )

    def time_polygon_query_loop(self, n_shapes, n_boxes):
        for poly in self.box_polygons:
            polygon_query(
                self.circles,
                polygon=poly,
                target_coordinate_system="global",
            )
