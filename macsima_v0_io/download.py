import os
from pathlib import Path
import subprocess

# Set to True to download only 3 sample images, False to download full dataset
DOWNLOAD_SAMPLE_ONLY = True

data_dir = Path(__file__).resolve().parent / "data"
os.makedirs(data_dir, exist_ok=True)

# Zenodo record for OMAP-10 (macsima v0 format)
ZENODO_RECORD = "7875938"
ARCHIVE_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD}/files-archive"
ARCHIVE_FILENAME = "OMAP10.zip"

# Sample files to download when DOWNLOAD_SAMPLE_ONLY is True
SAMPLE_FILES = [
    "DAPI_01.tif",
    "CD3_REAL1097.tif",
    "CD20_REA1087.tif",
]

if DOWNLOAD_SAMPLE_ONLY:
    # Download only 3 sample images
    print("NOTE: Downloading only a subset of the data (3 sample images).")
    print("The full OMAP-10 dataset contains 31 images (~6.1 GB).")
    print("Set DOWNLOAD_SAMPLE_ONLY = False to download the complete dataset.\n")
    for filename in SAMPLE_FILES:
        file_url = f"https://zenodo.org/records/{ZENODO_RECORD}/files/{filename}?download=1"
        output_path = data_dir / filename
        print(f"Downloading {filename}...")
        command = f'curl -L -o "{output_path}" "{file_url}"'
        subprocess.run(command, shell=True, check=True)
    print(f"Downloaded {len(SAMPLE_FILES)} sample files to {data_dir}")
else:
    # Download the full dataset archive
    command = f'curl -o "{data_dir / ARCHIVE_FILENAME}" "{ARCHIVE_URL}"'
    subprocess.run(command, shell=True, check=True)
    print(data_dir)
    # Unzip the data
    command = f'unzip -o "{data_dir / ARCHIVE_FILENAME}" -d "{data_dir}"'
    subprocess.run(command, shell=True, check=True)

