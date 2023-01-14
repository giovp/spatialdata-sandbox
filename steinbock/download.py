import os
from pathlib import Path
from tqdm import tqdm

urls = [
    "https://zenodo.org/record/7412972/files/cells.h5ad",
    "https://zenodo.org/record/7412972/files/masks_ilastik.zip",
    "https://zenodo.org/record/7412972/files/masks_deepcell.zip",
    "https://zenodo.org/record/7412972/files/img.zip",
    "https://zenodo.org/record/7412972/files/images.csv",
    "https://zenodo.org/record/7412972/files/ome.zip",
    "https://zenodo.org/record/7412972/files/panel.csv",
]

# os.makedirs("data", exist_ok=True)
# os.makedirs("data/steinbock", exist_ok=True)
# for url in tqdm(urls, desc="downloading"):
#     command = f"curl {url} --output 'data/steinbock/{Path(url).name}'"
#     os.system(command)

os.chdir("data/steinbock")
os.system("unzip masks_ilastik.zip")
os.system("unzip masks_deepcell.zip")
os.system("unzip ome.zip")

# os.chdir("data/xenium")
# os.system("unzip Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip")
