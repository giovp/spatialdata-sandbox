##
from spatialdata_io import visium_hd
import spatialdata as sd

##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("visium_hd_3.0.0_crc2"):
    path /= "visium_hd_3.0.0_crc2"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"

##
sdata = visium_hd(
    path_read,
    bin_size=8
)

sdata = sdata.filter_by_coordinate_system("downscaled_hires")

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