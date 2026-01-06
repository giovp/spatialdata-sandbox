import os
from pathlib import Path
import subprocess

data_dir = Path(__file__).resolve().parent / "data"
os.makedirs(data_dir, exist_ok=True)


url = "https://zenodo.org/api/records/14008816/files-archive"
filename = "OMAP23.zip"


# Download the data
command = f"curl -o {data_dir / filename} {url}"
subprocess.run(command, shell=True, check=True)

# Unzip the data
command = f"unzip -o {data_dir / filename} -d {data_dir}"
subprocess.run(command, shell=True, check=True)

