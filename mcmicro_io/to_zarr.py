from spatialdata_io import mcmicro
import spatialdata as sd
from pathlib import Path
import shutil

path = Path().resolve()
if not str(path).endswith("mcmicro_io"):
    path /= "mcmicro_io"
    assert path.exists()
path_read = path / "data" / "mcmicro"
path_write = path / "data.zarr"
##
print("parsing the data... ", end="")
sdata = mcmicro(
    path=path_read,
    dataset_id="exemplar-001"
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
