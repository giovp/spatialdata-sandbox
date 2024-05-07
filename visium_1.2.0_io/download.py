#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

##

# from https://www.10xgenomics.com/datasets/human-breast-cancer-targeted-immunology-panel-1-standard-1-2-0
# Human Breast Cancer: Targeted, Immunology Panel

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_BreastCancer_Immunology/Targeted_Visium_Human_BreastCancer_Immunology_image.tif",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_BreastCancer_Immunology/Targeted_Visium_Human_BreastCancer_Immunology_molecule_info.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_BreastCancer_Immunology/Targeted_Visium_Human_BreastCancer_Immunology_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Targeted_Visium_Human_BreastCancer_Immunology/Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz",
]

##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

subprocess.run("tar -xzf data/Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz -C data", shell=True, check=True)
subprocess.run("rm data/Targeted_Visium_Human_BreastCancer_Immunology_spatial.tar.gz", shell=True, check=True)
