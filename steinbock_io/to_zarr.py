from spatialdata_io import steinbock
import spatialdata as sd
from pathlib import Path
import shutil

path = Path().resolve()
if not str(path).endswith("steinbock_io"):
    path /= "steinbock_io"
    assert path.exists()
path_read = path / "data" / "steinbock"
path_write = path / "data.zarr"
##
print("parsing the data... ", end="")
sdata = steinbock(
    path=path_read,
)
print("done")
##
print("writing the data... ", end="")
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
##
print(f'view with "python -m napari_spatialdata view data.zarr"')
sdata = sd.SpatialData.read("./data.zarr/")
print("read")
print(sdata)
##
