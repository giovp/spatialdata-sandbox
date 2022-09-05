##
import json
import numpy as np
import scanpy as sc
import anndata as ad
import imageio
import shutil
from pathlib import Path
import spatialdata as sd

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
adata = sc.read_h5ad(path_read / "single_molecule.h5ad")
single_molecule = ad.AnnData(shape=(len(adata), 0))
single_molecule.obsm['spatial'] = adata.X
single_molecule.obsm['spatial_type'] = adata.obsm['cell_type']

j = json.load(open(path_read / "image_transform.json", "r"))
image_translation = np.array([j["translation_x"], j["translation_y"]])
image_scale_factors = np.array([j["scale_factor_x"], j["scale_factor_y"]])

expression = cells.copy()
del expression.obsm["region_radius"]
del expression.obsm["spatial"]
expression.obs['regions_key'] = 'cells'
expression.obs['instance_key'] = 'cell_id'
expression.obs['cell_id'] = np.arange(len(expression))

regions = ad.AnnData(shape=(len(cells.obsm['spatial']), 0))
regions.obsm['region_radius'] = cells.obsm['region_radius']
regions.obsm['spatial'] = cells.obsm['spatial']
regions.obs['cell_id'] = np.arange(len(regions))

##
transform = sd.Transform(translation=image_translation, scale_factors=image_scale_factors)

sdata = sd.SpatialData(
    table=expression,
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
print('view with "python -m napari_spatialdata view spatialdata-sandbox/merfish/data.zarr"')
##

