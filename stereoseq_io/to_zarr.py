##
from pathlib import Path
from spatialdata_io import stereoseq
import shutil
import spatialdata as sd

##
path = Path().resolve()
if not str(path).endswith("stereoseq_io"):
    path /= "stereoseq_io"
    assert path.exists()

# you can use symlinks to make the data available
# C://STT1_image_alignment//pipeline output//result"
dataset = Path("data/STT1_image_alignment/pipeline output/result")

path_read = path / dataset
sdata = stereoseq(path_read)

##
path_write = path / "data.zarr"
print("writing the data... ", end="")
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")

##
print(f'view with "python -m napari_spatialdata view data.zarr"')
sdata = sd.SpatialData.read(path_write)
print("read")
print(sdata)

##
