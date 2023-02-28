##
import os

os.environ["USE_PYGEOS"] = "0"
import json
import numpy as np
import scanpy as sc
import shutil
from pathlib import Path
import spatialdata as sd
import imageio.v3 as iio
import pandas as pd
from spatialdata._core.transformations import Scale, Translation, Sequence

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

translation = Translation(translation=image_translation, axes=("y", "x"))
scale = Scale(scale=image_scale_factors, axes=("y", "x"))
composed = Sequence([scale, translation])

img = iio.imread(path_read / "image.png")
img = np.expand_dims(img, axis=0)
img = sd.Image2DModel.parse(img, dims=("c", "y", "x"), transformations={'global': composed})
##
annotations = pd.DataFrame({"cell_type": pd.Categorical(adata.obsm["cell_type"])})
single_molecule = sd.PointsModel.parse(adata.X, annotation=annotations, feature_key="cell_type")

expression = cells.copy()
del expression.obsm["region_radius"]
del expression.obsm["spatial"]
expression.obs["cell_id"] = np.arange(len(cells))
expression = sd.TableModel.parse(
    adata=expression,
    region="/shapes/cells",
    instance_key="cell_id",
)
xy = cells.obsm["spatial"]
regions = sd.ShapesModel.parse(
    xy,
    geometry=0,
    radius=cells.obsm["region_radius"],
    index=expression.obs['cell_id'].copy()
)

polygons = sd.ShapesModel.parse(
    path_read / "anatomical.geojson"
)

##
sdata = sd.SpatialData(
    table=expression,
    shapes={"cells": regions, "anatomical": polygons},
    points={"single_molecule": single_molecule},
    images={"rasterized": img},
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
