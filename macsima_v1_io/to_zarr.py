from spatialdata_io import macsima
import spatialdata as sd

import shutil
from pathlib import Path

path_read = Path(__file__).resolve().parent / "data"
path_write = Path(__file__).resolve().parent / "data.zarr"

print("Parsing the data...")
sdata = macsima(path_read)
print("done")

print("Writing the data...")
if path_write.exists():
    shutil.rmtree(path_write)

sdata.write(path_write)
print("done")

sdata = sd.SpatialData.read("./data.zarr/")
print(sdata)