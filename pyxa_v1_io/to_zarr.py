# /// script
# requires-python = ">=3.12"
# dependencies = [
#     # the pyxa reader is pending review in scverse/spatialdata-io
#     "spatialdata-io @ git+https://github.com/ckmah/spatialdata-io.git@pyxa-reader",
# ]
# ///
##
from spatialdata_io.experimental import pyxa
import spatialdata as sd

##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("pyxa_v1_io"):
    path /= "pyxa_v1_io"
    assert path.exists()

path_read = path / "data"
path_write = path / "data.zarr"

##
print("parsing the data... ", end="")
sdata = pyxa(path_read, image_path=path_read / "mosaic_3d.ome.zarr")
print("done")

##
print("writing the data... ", end="")
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")

##
sdata = sd.SpatialData.read(path_write)
print(sdata)
