##
from spatialdata_io import visium
import spatialdata as sd
##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("visium2"):
    path /= "visium2"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"
##
sdata = visium(path_read)
##
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
##
print(f'view with "python -m spatialdata view data.zarr"')
sdata = sd.SpatialData.read(path_write)
print(sdata)
print("read")
