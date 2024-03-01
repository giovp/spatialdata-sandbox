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
sdataST8059048 = visium(
    path_read / "mouse_brain_visium_wo_cloupe_data/rawdata/ST8059048",
    dataset_id="ST8059048",
)
sdataST8059050 = visium(
    path_read / "mouse_brain_visium_wo_cloupe_data/rawdata/ST8059050",
    dataset_id="ST8059050",
)
# Each SpatialData object has 3 coordinate systems: 'global', 'downscaled_hires' and 'downscaled_lowres'.
# Visium datasets generally are available with images (lowres and hires) and sometimes also a "full resolution" image,
# which has an even greater resolution than the hires image. This dataset doesn't contain a full resolution image, so
# let's drop the coordinate system "global", which would otherwise contain it. Let's also drop the coordinate system
# "downscaled_lowres", which is not needed.
# NOTE: in future version of the Visium reader only one coordinate system will be used, containing all the images. We
# are keeping the three coordinate systems for legacy reasons for the time being.
for sdata in [sdataST8059048, sdataST8059050]:
    for el in sdata._gen_spatial_element_values():
        for cs_name in ["global", "downscaled_lowres"]:
            if cs_name in sd.transformations.get_transformation(el, get_all=True):
                sd.transformations.remove_transformation(el, cs_name)

# we want to concatenate the two datasets, but they have the same coordinate system names; therefore let's rename
# the each coordinate system to match the dataset id 
sdataST8059048.rename_coordinate_systems(
    {"downscaled_hires": "ST8059048"}
)
sdataST8059050.rename_coordinate_systems(
    {"downscaled_hires": "ST8059050"}
)
sdata = sd.concatenate([sdataST8059048, sdataST8059050], concatenate_tables=True)
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
