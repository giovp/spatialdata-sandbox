## loading multiple samples visium data from disk (SpaceRanger), concatenating and saving them to .zarr

import spatialdata as sd
from spatialdata_io import read_visium

samples = ["152806", "152807", "152810", "152811"]
sdatas = []

for sample in samples:
    sdata = read_visium(path=sample, coordinate_system_name=sample)
    sdatas.append(sdata)

sdata = sd.SpatialData.concatenate(sdatas, merge_tables=True)
sdata.write("data.zarr")

## calling the SpatialData constructor with some transformations on it
# REAL CODE!
# load arbitarily-stored data (in this case: .h5ad, .geojson and .png files)
import numpy as np
import scanpy as sc
import spatialdata as sd
import imageio.v3 as iio
import xarray as xr
from spatialdata_io import (
    points_anndata_from_coordinates,
    polygons_anndata_from_geojson,
    circles_anndata_from_coordinates,
    table_update_anndata,
)

cells = sc.read_h5ad("cells.h5ad")
img = xr.DataArray(iio.imread("image.png"), dims=("y", "x"))
adata = sc.read_h5ad("single_molecule.h5ad")

# construct spatial elements
# points
single_molecule = points_anndata_from_coordinates(coordinates=adata.X, points_types=adata.obsm["cell_type"])

# circles
xy = cells.obsm["spatial"]
regions = circles_anndata_from_coordinates(
    coordinates=xy,
    radii=cells.obsm["region_radius"],
    instance_key="cell_id",
    instance_values=np.arange(len(xy)),
)

# polygons
adata_polygons = polygons_anndata_from_geojson("anatomical.geojson")

# transformations
translation = sd.Translation(translation=np.array([123., 10.]))
scale = sd.Scale(scale=np.array([0.5, 0.5]))
composed = sd.Sequence([scale, translation])

# annotation table (cell expression matrix)
expression = cells.copy()
del expression.obsm["region_radius"]
del expression.obsm["spatial"]
table_update_anndata(
    adata=expression,
    regions="/points/cells",
    regions_key="regions_id",
    instance_key="cell_id",
    regions_values="/points/cells",
    instance_values=np.arange(len(cells)),
)
# constructor
sdata = sd.SpatialData(
    table=expression,
    points={"cells": regions, "single_molecule": single_molecule},
    images={"rasterized": img},
    polygons={"anatomical": adata_polygons},
    transformations={
        ("/images/rasterized", "global"): composed,
    },
)
# limitations
# - no support for transformations from intrinsitc element to intrinstic element. Do we want them?
# - shall we allow for explicly specifying the coordinate systems?
# - should transformations be defined outside the spatialdata data constructor?

# PSEUDOCODE!
# annotation table doesn't change
# points, circles, polygons constructions don't change
# images
img = xr.DataArray(iio.imread("image.png"), dims=("y", "x"))

cs = sd.CoordinateSystem(name="global", axes=["x", "y", 'ChaNnEls'])

image = xr.DataArray(...)
translation = sd.Translation(translation=np.array([123., 10.]))
scale = sd.Scale(scale=np.array([0.5, 0.5]))
composed = sd.Sequence([scale, translation])
composed.input = '/images/my_image'
composed.output = cs

image_add_transformation(xarray, transformation=composed)
image.transformations['my_cs']

image.commit_transformation(target_cs='my_cs')
tr = image.transformations['my_cs']
tr.apply(image.data)


# constructor
sdata = sd.SpatialData(
    table=expression,
    points={"cells": regions, "single_molecule": single_molecule},
    images={"rasterized": img},
    polygons={"anatomical": adata_polygons},
    transformations={
        ("extrinsic0", "extrinsic1"): composed
    }
)
sdata.to_zarr(...)

new_image = xr.DataArrat()
image_add_transformation(new_image, transfomation=sdata.images['my_image'])
sdata['new_image'] = new_image

## aggregation with multiple types of elements
# MORE DISCUSSION NEEDED
sdata = from_zarr('..')

img = sdata.images['my_image']
visium_sposts = sdata.circles['visium']
# source = sdata.points['transcripts']
# target = sdata.polygons['anatomical']

# RETURN TYPE
# table = spatialdata.aggregate(source=img, regions=visium_spots, method: Union[List[str], List[Callable]] =[np.mean, np.std])  #
# method =
# 'mean'
sdata.aggregate(source='/img', regions='/sample0/circles/visium')


# not supported
# img.aggregate_by(regions=visium_spots)
# visium_sposts.aggregate(source=img)

## loading specific samples from a multi-sample dataset (subsetting by coordinate system)
# reader with a flag to parse and infer the hierarcy from the ngff file into coordinate systems
"""
SpatialData object with:
├── images
│     ├── '/images/point16': DataArray (3, 1024, 1024), with axes: c, y, x
│     ├── '/images/point23': DataArray (3, 1024, 1024), with axes: c, y, x
│     └── 'point8': DataArray (3, 1024, 1024), with axes: c, y, x
├── labels
│     ├── '/labels/point16': DataArray (1024, 1024), with axes: y, x
│     ├── 'point23': DataArray (1024, 1024), with axes: y, x
│     └── 'point8': DataArray (1024, 1024), with axes: y, x
├── polygons
│     └── 'Shapes_point16_1': AnnData with obs.spatial describing 2 polygons, with axes x, y
└── table
      └── 'AnnData object with n_obs × n_vars = 3309 × 36
    obs: 'row_num', 'point', 'cell_id', 'X1', 'center_rowcoord', 'center_colcoord', 'cell_size', 'category', 'donor', 'Cluster', 'batch', 'library_id'
    uns: 'mapping_info'
    obsm: 'X_scanorama', 'X_umap', 'spatial'': AnnData (3309, 36)
with coordinate systems:
▸ point16
    with axes: c, y, x
    with elements: /images/point16, /labels/point16, /polygons/Shapes_point16_1
▸ point23
    with axes: c, y, x
    with elements: /images/point23, /labels/point23
▸ point8
    with axes: c, y, x
    with elements: /images/point8, /labels/point8
"""
sdata0 = sdata.query.coordinate_system(samples[0], filter_rows=False)
sdata1 = sdata.query.bounding_box()
sdata1 = sdata.query.polygon('/polygons/annotations')
sdata1 = sdata.query.table