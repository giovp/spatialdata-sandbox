##
import numpy as np
import anndata as ad
import shutil
from pathlib import Path
import spatialdata as sd
from spatialdata.transformations.transformations import Scale, Identity
import re
import os
from tqdm import tqdm
from spatialdata_io import visium

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("visium"):
    path /= "visium"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"

##
sdata = visium(path_read / 'mouse_brain_visium_wo_cloupe_data/rawdata/ST8059048', dataset_id='ST8059048')
print(sdata)

#
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
print(f'view with "python -m napari_spatialdata view data.zarr"')
sdata = sd.SpatialData.read(path_write)
print(sdata)
print("read")
