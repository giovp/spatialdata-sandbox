import os
import subprocess
from pathlib import Path
import shutil

# URL of the dataset
url = "https://s3.embl.de/spatialdata/raw_data/mouseLiver.zarr.zip"

# Create the data directory if it doesn't exist
data_dir = Path("data")
data_dir.mkdir(parents=True, exist_ok=True)

# Download the dataset
zip_path = data_dir / "mouseLiver.zarr.zip"
subprocess.run(["curl", "-o", zip_path, url], check=True)

print("Download complete.")
