##
from spatialdata_io import xenium
import spatialdata as sd

##
from pathlib import Path
import shutil

##
path = Path().resolve()
# luca's workaround for pycharm
if not str(path).endswith("xenium_rep1_io"):
    path /= "xenium_rep1_io"
    assert path.exists()
# new path
# path_read = path / "data/xenium"
# old path
path_read = path / "data/xenium/outs"
path_write = path / "data.zarr"
##
print("parsing the data... ", end="")
sdata = xenium(
    path=str(path_read),
    n_jobs=8,
    cells_as_shapes=True,
    # morphology_mip=False,
    # morphology_focus=False,
    # nucleus_boundaries=False,
    # cell_boundaries=False,
    # points=False
)
print("done")
##
print("writing the data... ", end="")
if path_write.exists():
    shutil.rmtree(path_write)
sdata.write(path_write)
print("done")
##
print(f'view with "python -m spatialdata view data.zarr"')
sdata = sd.SpatialData.read("./data.zarr/")
print("read")
print(sdata)
##
