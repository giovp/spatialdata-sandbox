##
import numpy as np
import scanpy as sc
import os
import imageio
from pathlib import Path
import spatialdata as sd
import ngff_tables_prototype as viz

##
print(f"os.getcwd() = {os.getcwd()}")
data_dir = Path().resolve() / "merfish_mouse_brain" / "processed"
assert data_dir.exists()

##
cells = sc.read_h5ad(data_dir / "cells.h5ad")
img = np.asarray(imageio.imread(data_dir / "image.png"))
sm = sc.read_h5ad(data_dir / "single_molecule.h5ad")

expression = cells.copy()
del expression.obsm['region_radius']
del expression.obsm['spatial']

regions = cells.copy()
del regions.X

sdata = sd.SpatialData(adata=expression, regions={'cells': regions}, images={'rasterized': img}, points=sm)

##
sdata.to_zarr(data_dir / "spatialdata.zarr")