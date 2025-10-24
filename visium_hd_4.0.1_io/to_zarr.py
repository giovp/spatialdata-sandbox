##
from spatialdata_io import visium_hd
import spatialdata as sd

##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("visium_hd_4.0.1_io"):
    path /= "visium_hd_4.0.1_io"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"

##
sdata = visium_hd(
    path_read,
    load_segmentations_only=True,
    # needed for centroids but otherwise don't have to load nuclei
    # Will increase load time as matrix is calculated from the 2um binned data
    load_nucleus_segmentations=True,
)

##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")

##
print(f'view with "python -m napari_spatialdata view data.zarr"')
sdata = sd.SpatialData.read(path_write)
print(sdata)
print("read")

# from napari_spatialdata import Interactive
# Interactive(sdata)
