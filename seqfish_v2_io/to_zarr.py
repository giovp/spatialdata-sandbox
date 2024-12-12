from spatialdata_io import seqfish
import spatialdata as sd
from pathlib import Path
import shutil

path = Path().resolve()
if not str(path).endswith("seqfish_v2_io"):
    path /= "seqfish_v2_io"
    assert path.exists()
path_read = path / "data"
path_write = path / "data.zarr"
##
print("parsing the data... ", end="")
sdata = seqfish(
    path=path_read,
    # load_images=True,
    # load_labels=True,
    # load_points=True,
    # rois=[1],
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

from napari_spatialdata import Interactive
Interactive(sdata)