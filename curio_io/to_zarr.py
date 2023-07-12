##
from pathlib import Path
from spatialdata_io import curio
import shutil
import spatialdata as sd

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("curio_io"):
    path /= "curio_io"
    assert path.exists()

##
# Mouse_hippocampus_output appears to be incomplete, excluding it for now
# datasets = ["Mouse_liver_output", "Mouse_spleen_output"]
datasets = ["Mouse_liver_output"]

sdatas = []
for dataset in datasets:
    path_read = path / "data" / dataset
    sdata = curio(path_read)
    # from napari_spatialdata import Interactive
    # Interactive(sdata)
    sdatas.append(sdata)
##
# we will write to Zarr only the first dataset
path_write = path / "data.zarr"
sdata = sdatas[0]

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
