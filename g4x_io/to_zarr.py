##
from spatialdata_io import g4x
import spatialdata as sd

##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("g4x_io"):
    path /= "g4x_io"
    assert path.exists()

path_read = path / "data"
path_write = path / "data.zarr"

##
print("parsing the data... ", end="")
sdata = g4x(
    input_path=str(path_read),
    output_path=str(path_write),
    include_he=True,
    include_segmentation=True,
    include_protein=True,
    include_transcripts=True,
    include_tables=True,
    mode="append",
)
print("done")

##
sdata = sd.SpatialData.read("./data.zarr/")
print(sdata)
##
