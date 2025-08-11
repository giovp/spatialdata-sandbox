import subprocess
from pathlib import Path
import shutil
import spatialdata as sd

data_dir = Path("data")
zip_path = data_dir / "mouseLiver.zarr.zip"

# Unzip the dataset into the extraction directory
subprocess.run(["unzip", "-o", zip_path, "-d", '.'], check=True)

# reading and rewriting ensures we are writing according to the latest format
src_path = Path('mouseLiver.zarr')
des_path = Path('data.zarr')

sdata = sd.read_zarr(src_path)
print(sdata)

if src_path.exists():
    shutil.rmtree(des_path)

sdata.write(des_path)

# cleanup
if src_path.exists():
    shutil.rmtree(src_path)

print("Write to Zarr complete.")
