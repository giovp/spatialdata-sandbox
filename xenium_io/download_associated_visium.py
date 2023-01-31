##
import os
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
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_analysis.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer"
    "/CytAssist_FFPE_Human_Breast_Cancer_molecule_info.h5",
]

##
os.makedirs("data", exist_ok=True)
os.makedirs("data/visium", exist_ok=True)
os.chdir('data/visium')
for url in tqdm(urls, desc="downloading"):
    command = f"curl -O {url} --output {'data/visium/' + Path(url).name}"
    if not url.endswith(".tar.gz"):
        os.system(command)

##
os.system('tar -xvf CytAssist_FFPE_Human_Breast_Cancer_analysis.tar.gz')
os.system('tar -xvf CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz')
os.system('mv CytAssist_FFPE_Human_Breast_Cancer_tissue_image.tif CytAssist_FFPE_Human_Breast_Cancer_image.tif')
##
# from spatialdata_io import read_visium
# sdata = read_visium('.')