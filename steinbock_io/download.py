import os
from pathlib import Path
from tqdm import tqdm
import subprocess

urls = [
    "https://zenodo.org/record/7412972/files/cells.h5ad",
    "https://zenodo.org/record/7412972/files/masks_ilastik.zip",
    "https://zenodo.org/record/7412972/files/masks_deepcell.zip",
    "https://zenodo.org/record/7412972/files/img.zip",
    "https://zenodo.org/record/7412972/files/images.csv",
    "https://zenodo.org/record/7412972/files/ome.zip",
    "https://zenodo.org/record/7412972/files/panel.csv",
]

os.makedirs("data", exist_ok=True)
os.makedirs("data/steinbock", exist_ok=True)
for url in tqdm(urls, desc="downloading"):
    command = f"curl {url} --output 'data/steinbock/{Path(url).name}'"
    subprocess.run(command, shell=True, check=True)

os.chdir("data/steinbock")
subprocess.run("unzip -f -o masks_ilastik.zip", shell=True, check=True)
subprocess.run("unzip -f -o masks_deepcell.zip", shell=True, check=True)
subprocess.run("unzip -f -o ome.zip", shell=True, check=True)
