import concurrent.futures
import requests
import zipfile
import os
from pathlib import Path

SAVE = Path().resolve() / "data" 
SAVE.mkdir(parents=True, exist_ok=True)
assert SAVE.exists()
raw = SAVE / "raw"
raw.mkdir(parents=True, exist_ok=True)
processed = SAVE / "processed"
processed.mkdir(parents=True, exist_ok=True)

def download_file(url, output_dir):
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        file_name = output_dir / Path(url).name

        with open(file_name, "wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    file.write(chunk)

        print(f"Downloaded {url} to {file_name}")

    except Exception as e:
        print(f"Error downloading {url}: {e}")

def extract_zip(zip_path, extract_path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_file_list = zip_ref.namelist()
        zip_ref.extractall(extract_path, members=[f for f in zip_file_list if not f.startswith("__MACOSX/")])

if __name__ == "__main__":
    urls = [
        "http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom",
        "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/data/mRNA_coords_raw_counting.hdf5",
        "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/data/polyT_seg.pkl",
        "https://s3.embl.de/spatialdata/raw_data/im_nuc_small.pickle"
    ]

    p_urls = [
        "https://s3.embl.de/spatialdata/raw_data/ssamdataset-osmFISH.zarr.zip"
    ]


    output_directory = raw
    p_output_directory = processed

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = []
        for url in urls:
            future = executor.submit(download_file, url, output_directory)
            futures.append(future)
        for url in p_urls:
            future = executor.submit(download_file, url, p_output_directory)
            futures.append(future)
        concurrent.futures.wait(futures)

    for p_url in p_urls:
        zip_filename = p_output_directory / Path(p_url).name
        extract_path = p_output_directory
        if zip_filename.exists() and zip_filename.suffix == ".zip":
            extract_zip(zip_filename, extract_path)

    