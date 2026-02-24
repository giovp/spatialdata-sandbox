import os
from pathlib import Path
import subprocess

# Set to True to download only 3 sample images, False to download full dataset
DOWNLOAD_SAMPLE_ONLY = True

data_dir = Path(__file__).resolve().parent / "data"
os.makedirs(data_dir, exist_ok=True)

# Zenodo record for OMAP-23 (macsima v1 format)
ZENODO_RECORD = "14008816"
ARCHIVE_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD}/files-archive"
ARCHIVE_FILENAME = "OMAP23.zip"

# Sample files to download when DOWNLOAD_SAMPLE_ONLY is True
# Most images are .ome.tiff, but there is a file with extension .tif. 
# Let's download this case to cover both cases.
# The single .png image is just as a cover image for Zenodo and not representative.
SAMPLE_FILES = [
    "OMAP23_DAPI_01.ome.tif",
    "OMAP23_CD3E.ome.tif",
    "OMAP23_FoxP3.tif",
]

if DOWNLOAD_SAMPLE_ONLY:
    # Download only 3 sample images
    print("NOTE: Downloading only a subset of the data (4 sample images).")
    print("The full OMAP-23 dataset contains 30 images (~5.6 GB).")
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
    # Unzip the data
    command = f'unzip -o "{data_dir / ARCHIVE_FILENAME}" -d "{data_dir}"'
    subprocess.run(command, shell=True, check=True)

