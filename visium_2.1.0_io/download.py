#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

##

# from https://www.10xgenomics.com/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard
# Visium CytAssist Gene and Protein Expression Library of Human Tonsil, H&E, 6.5 mm (FFPE)

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_tissue_image.tif",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_molecule_info.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_spatial.tar.gz",
]

##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

subprocess.run("tar -xzf data/CytAssist_FFPE_Protein_Expression_Human_Tonsil_spatial.tar.gz -C data", shell=True, check=True)
subprocess.run("rm data/CytAssist_FFPE_Protein_Expression_Human_Tonsil_spatial.tar.gz", shell=True, check=True)
