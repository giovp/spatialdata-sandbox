##
import numpy as np
from dask_image.ndinterp import affine_transform as affine_transform_dask
from skimage import data
from xarray import DataArray
import dask.array as da

image = data.astronaut()
image = DataArray(image, dims=["y", "x", "c"])

scale_factor = 10
matrix = np.array(
    [
        [1.0 / scale_factor, 0, 0, 0],
        [0, 1.0 / scale_factor, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
)
output_shape = (image.shape[0] * scale_factor, image.shape[1] * scale_factor, 3)
image_transformed = affine_transform_dask(
    image.data,
    matrix=matrix,
    output_shape=output_shape,
)
image_transformed = DataArray(image_transformed, dims=("y", "x", "c"))
image_transformed
##
import spatialdata as sd
from napari_spatialdata import Interactive

image_transformed0 = sd.models.Image2DModel.parse(
    image_transformed,
    dims=("y", "x", "c"),
    transformations={"global": sd.transformations.Scale([0.5, 1.1], axes=("x", "y"))},
)
image_transformed1 = sd.models.Image2DModel.parse(
    image_transformed.copy(),
    dims=("y", "x", "c"),
    transformations={"global": sd.transformations.Translation([500, -1000], axes=("x", "y"))},
)
image_transformed2 = sd.models.Image2DModel.parse(
    image_transformed.copy(),
    dims=("y", "x", "c"),
    transformations={
        "global": sd.transformations.Affine(
            [[0.25, 0.5, -500], [1, 0.125, -1000], [0, 0, 1]], input_axes=("x", "y"), output_axes=("x", "y")
        )
    },
)
sdata = sd.SpatialData(images={"im2": image_transformed2, "im0": image_transformed0, "im1": image_transformed1})
print(sdata)

# TODO: rotate only the part with determinant = 1
t = image_transformed2.attrs["transform"]["global"]
image_transformed3 = sd.transform(image_transformed2, t, maintain_positioning=True)
sdata["im3"] = image_transformed3

import pandas as pd
from geopandas import GeoDataFrame
from shapely import Polygon
from spatialdata._core.query.spatial_query import get_bounding_box_corners

min_coordinates_data_bb = np.array([0, 0])
max_coordinates_data_bb = np.array([image_transformed.shape[1], image_transformed.shape[0]])
data_bb = get_bounding_box_corners(
    axes=("x", "y"), min_coordinate=min_coordinates_data_bb, max_coordinate=max_coordinates_data_bb
)
data_pol = Polygon(data_bb)


def parse_polygon(polygons: list[Polygon]) -> GeoDataFrame:
    return sd.models.ShapesModel.parse(
        GeoDataFrame(geometry=polygons), transformations={"global": sd.transformations.Identity()}
    )


data_pol0 = sd.transform(parse_polygon([data_pol]), image_transformed0.attrs["transform"]["global"])
data_pol1 = sd.transform(parse_polygon([data_pol]), image_transformed1.attrs["transform"]["global"])
data_pol2 = sd.transform(parse_polygon([data_pol]), image_transformed2.attrs["transform"]["global"])

min_coordinates_query_bb = np.array([-850, -750])
max_coordinates_query_bb = np.array([1000, 1000])
query_bb = get_bounding_box_corners(
    axes=("x", "y"), min_coordinate=min_coordinates_query_bb, max_coordinate=max_coordinates_query_bb
)
query_pol = Polygon(query_bb)
bounding_boxes = pd.concat([data_pol0, data_pol1, data_pol2, parse_polygon([query_pol])])
##
TARGET_WIDTH = 1000
rasterized_sdata = sd.rasterize(
    sdata,
    axes=("x", "y"),
    min_coordinate=min_coordinates_query_bb,
    max_coordinate=max_coordinates_query_bb,
    target_coordinate_system="global",
    target_width=TARGET_WIDTH,
)

import rasterio.features
from spatialdata._core.operations.rasterize import _compute_target_dimensions

target_width, target_height, target_depth = _compute_target_dimensions(
    spatial_axes=("x", "y"),
    min_coordinate=min_coordinates_query_bb,
    max_coordinate=max_coordinates_query_bb,
    target_unit_to_pixels=None,
    target_width=TARGET_WIDTH,
    target_height=None,
    target_depth=None,
)
target_unit_to_pixels_x = target_width / (max_coordinates_query_bb - min_coordinates_query_bb)[0]
target_unit_to_pixels_y = target_height / (max_coordinates_query_bb - min_coordinates_query_bb)[1]

##
# translations in units
translation = sd.transformations.Translation(min_coordinates_query_bb, axes=("x", "y"))
# pixels to units
scale = sd.transformations.Scale([target_unit_to_pixels_x, target_unit_to_pixels_y], axes=("x", "y")).inverse()
# needs to be the transformation from pixels to units
sequence = sd.transformations.Sequence([scale, translation])
rasterio_affine = rasterio.transform.Affine(
    *sequence.to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))[:2, :].flatten().tolist()
)

data_pol0_raster = rasterio.features.rasterize(
    shapes=data_pol0.geometry, out_shape=(int(target_height), int(target_width)), transform=rasterio_affine
)
data_pol1_raster = rasterio.features.rasterize(
    shapes=data_pol1.geometry, out_shape=(int(target_height), int(target_width)), transform=rasterio_affine
)
data_pol2_raster = rasterio.features.rasterize(
    shapes=data_pol2.geometry, out_shape=(int(target_height), int(target_width)), transform=rasterio_affine
)
import matplotlib.pyplot as plt

# plt.figure()
# plt.imshow(data_pol0_raster)
# plt.show()


##
def stitch(target_data: np.array, source_image: np.array, mask: np.array) -> None:
    axes = plt.subplots(1, 4, figsize=(20, 5))[1].flatten()
    axes[0].imshow(target_data.transpose(1, 2, 0))
    axes[0].set_title("target_data before")

    axes[1].imshow(source_image.transpose(1, 2, 0))
    axes[1].set_title("source_image")

    axes[2].imshow(mask)
    axes[2].set_title("mask")

    y, x = np.where(mask)
    # y, x = da.where(mask)
    target_data[:, y, x] = source_image[:, y, x]

    axes[3].imshow(target_data.transpose(1, 2, 0))
    axes[3].imshow(source_image.transpose(1, 2, 0))
    plt.show()


# if we remove "compute()" we get a dask not implemented error since nd "sparse" indices are not supported
rasterized_image0 = rasterized_sdata["im0_rasterized_images"].data.compute()
rasterized_image1 = rasterized_sdata["im1_rasterized_images"].data.compute()
rasterized_image2 = rasterized_sdata["im2_rasterized_images"].data.compute()

target_data = np.zeros_like(rasterized_image0)
stitch(target_data=target_data, source_image=rasterized_image2, mask=data_pol2_raster.astype(bool))
stitch(target_data=target_data, source_image=rasterized_image0, mask=data_pol0_raster.astype(bool))
stitch(target_data=target_data, source_image=rasterized_image1, mask=data_pol1_raster.astype(bool))

plt.figure()
plt.imshow(target_data.transpose(1, 2, 0))
plt.show()
##


Interactive([rasterized_sdata, sd.SpatialData(shapes={"bounding boxes": bounding_boxes})])


Interactive([sdata, sd.SpatialData(shapes={"bounding boxes": bounding_boxes})])
