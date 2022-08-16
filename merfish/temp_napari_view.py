import shutil

import json
import numpy as np
import scanpy as sc
import os
import imageio
import shutil
from pathlib import Path
import spatialdata as sd
# import ngff_tables_prototype as viz
import napari_spatialdata as viz

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith('merfish'):
    path /= 'merfish'
path /= "data.zarr"
##
sdata = sd.SpatialData.read(path)
sdata
