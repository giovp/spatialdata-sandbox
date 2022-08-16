##
import shutil

import json
import numpy as np
import scanpy as sc
import os
import imageio
import shutil
from pathlib import Path
import spatialdata as sd
import ngff_tables_prototype as viz

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith('merfish'):
    path /= 'merfish'
##
path_read = path / "data" / "processed"
path_write = path / "data.zarr"

##
cells = sc.read_h5ad(path_read / "cells.h5ad")
img = np.asarray(imageio.imread(path_read / "image.png"))
single_molecule = sc.read_h5ad(path_read / "single_molecule.h5ad")
j = json.load(open(path_read / "image_transform.json", "r"))
image_translation = np.array([j["translation_x"], j["translation_y"]])
image_scale_factors = np.array([j["scale_factor_x"], j["scale_factor_y"]])

expression = cells.copy()
del expression.obsm["region_radius"]
del expression.obsm["spatial"]

regions = cells.copy()
del regions.X

##
transform = sd.Transform(translation=image_translation, scale_factors=image_scale_factors)

sdata = sd.SpatialData(
    tables={'cells': expression},
    points={"cells": regions, 'single_molecule': single_molecule},
    images={"rasterized": img},
    images_transform={"rasterized": transform},
)
print(sdata)

##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print('done')
##

