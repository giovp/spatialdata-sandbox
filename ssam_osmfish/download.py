import concurrent.futures
import requests
from pathlib import Path
import ssam
import pandas as pd

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
    

def download_ssam(spot):
    spots = pd.read_csv(spot, usecols=["xc", "yc", "gene", "cell"]).rename(columns={"xc": "x", "yc": "y"}).set_index('gene')
    ds = ssam.SSAMDataset("data/processed/ssam-osmfish")
    analysis = ssam.SSAMAnalysis(ds, ncores=40, verbose=True)

if __name__ == "__main__":
    urls = [
        "https://s3.amazonaws.com/starfish.data.spacetx/spacetx-website/data/smFISH_Allen/s3_cell_by_gene.csv",
        "https://s3.amazonaws.com/starfish.data.spacetx/spacetx-website/data/smFISH_Allen/s3_mapped_cell_table.csv",
        "https://s3.amazonaws.com/starfish.data.spacetx/spacetx-website/data/smFISH_Allen/s3_spot_table.csv",
    ]

    output_directory = raw

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        for url in urls:
            executor.submit(download_file, url, output_directory)