#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

##

# from https://www.10xgenomics.com/datasets/fresh-frozen-visium-on-cytassist-human-breast-cancer-probe-based-whole-transcriptome-profiling-2-standard
# Fresh Frozen Visium on CytAssist: Human Breast Cancer, Probe-Based Whole Transcriptome Profiling

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_image.tif",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_tissue_image.tif",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_molecule_info.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_Fresh_Frozen_Human_Breast_Cancer/CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz",
]

##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

subprocess.run("tar -xzf data/CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz -C data", shell=True, check=True)
subprocess.run("rm data/CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz", shell=True, check=True)
