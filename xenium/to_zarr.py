##
from spatialdata_io import convert_xenium_to_ngff
import spatialdata as sd
##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("xenium"):
    path /= "xenium"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"
##
# convert_xenium_to_ngff(path_read, path_write)
print("done")
##
print(f'view with "python -m spatialdata view data.zarr"')
print("read")
sdata = sd.SpatialData.read("./data.zarr/")
print(sdata)
