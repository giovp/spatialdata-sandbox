#!/usr/bin/env python3

from pathlib import Path
import spatialdata as sd

# Dataset is already in SpatialData format. Still, read it and convert it to the latest version
path_read = Path(__file__).parent / "data/data.zarr"
path_write = Path(__file__).parent / "data.zarr"

assert path_read.exists()
sdata = sd.SpatialData.read(path_read)

# Test reading
print(sdata)

sdata.write(path_write, overwrite=True)
print(f'view with "python -m napari_spatialdata view data.zarr"')