##
import json
import numpy as np
import scanpy as sc
import shutil
from pathlib import Path
import spatialdata as sd
import imageio.v3 as iio

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("merfish"):
    path /= "merfish"
    assert path.exists()
##
path_read = path / "data" / "processed"
path_write = path / "data.zarr"

##
cells = sc.read_h5ad(path_read / "cells.h5ad")
adata = sc.read_h5ad(path_read / "single_molecule.h5ad")
##
j = json.load(open(path_read / "image_transform.json", "r"))
image_translation = np.array([j["translation_y"], j["translation_x"]])
image_scale_factors = np.array([j["scale_factor_y"], j["scale_factor_x"]])

translation = sd.Translation(translation=image_translation)
scale = sd.Scale(scale=image_scale_factors)
composed = sd.Sequence([scale, translation])

img = iio.imread(path_read / "image.png")
img = np.expand_dims(img, axis=0)
img = sd.Image2DModel.parse(img, dims=("c", "y", "x"), transform=composed)
##
single_molecule = sd.PointsModel.parse(coords=adata.X, points_assignment=adata.obsm["cell_type"])

expression = cells.copy()
del expression.obsm["region_radius"]
del expression.obsm["spatial"]
expression = sd.TableModel.parse(
    adata=expression,
    region="/shapes/cells",
    instance_key="cell_id",
    instance_values=np.arange(len(cells)),
)
xy = cells.obsm["spatial"]
regions = sd.ShapesModel.parse(
    coords=xy,
    shape_type='Circle',
    shape_size=np.mean(cells.obsm["region_radius"]).item(),
    instance_key="cell_id",
    instance_values=np.arange(len(xy)),
)

with open(path_read / "anatomical.geojson") as infile:
    geojson_string = infile.read()
adata_polygons = sd.PolygonsModel.parse(geojson_string, instance_key="region_id")

##
sdata = sd.SpatialData(
    table=expression,
    shapes={"cells": regions},
    points={"single_molecule": single_molecule},
    images={"rasterized": img},
    polygons={"anatomical": adata_polygons},
)
print(sdata)
##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
print(f'view with "python -m napari_spatialdata view data.zarr"')
##
sdata = sd.SpatialData.read(path_write)
print(sdata)
print("read")
