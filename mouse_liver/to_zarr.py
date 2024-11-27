import os
import subprocess
from pathlib import Path
import shutil

data_dir = Path("data")
zip_path = data_dir / "mouseLiver.zarr.zip"

# Unzip the dataset into the extraction directory
subprocess.run(["unzip", "-o", zip_path, "-d", '.'], check=True)

shutil.move("mouseLiver.zarr", "data.zarr")

print("Extraction complete.")
