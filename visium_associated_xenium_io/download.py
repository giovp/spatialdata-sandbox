##
import os
import subprocess
import shutil
from pathlib import Path
from tqdm import tqdm

# the data that can be downloaded from the 10x website changed, this script downloads the newest version of
# the data

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_probe_set.csv",
    # this image has bad quality and it is not aligned. We need the tissue_image.tif
    # "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    # "/CytAssist_FFPE_Human_Breast_Cancer_image.tif",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_tissue_image.tif",
    # "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    # "/CytAssist_FFPE_Human_Breast_Cancer_analysis.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_molecule_info.h5",
]

##
os.makedirs("data", exist_ok=True)
for url in tqdm(urls, desc="downloading"):
    command = f"curl -o {'data/' + Path(url).name} {url}"
    subprocess.run(command, shell=True, check=True)

##
# subprocess.run("tar -xvf data/CytAssist_FFPE_Human_Breast_Cancer_analysis.tar.gz", shell=True, check=True)
subprocess.run("tar -xvf data/CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz -C data/", shell=True, check=True)
