#!/usr/bin/env python3

from pathlib import Path
import spatialdata as sd

# Dataset is already in SpatialData format.
path_read = Path(__file__).parent / "data.zarr"
assert path_read.exists()

print(f'view with "python -m napari_spatialdata view data.zarr"')

# Test reading
sdata = sd.SpatialData.read("./data.zarr")
print(sdata)
