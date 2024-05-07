#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

##

# from https://www.10xgenomics.com/datasets/visium-cytassist-gene-expression-libraries-of-post-xenium-human-colon-cancer-ffpe-using-the-human-whole-transcriptome-probe-set-2-standard
# Spatial Gene Expression Dataset by Space Ranger 2.1.0

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_tissue_image.btf",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_molecule_info.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_spatial.tar.gz",
]

##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

subprocess.run("tar -xzf data/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_spatial.tar.gz -C data", shell=True, check=True)
subprocess.run("rm data/CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_spatial.tar.gz", shell=True, check=True)
