#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

##

# from https://www.10xgenomics.com/datasets/gene-and-protein-expression-library-of-human-glioblastoma-cytassist-ffpe-2-standard
# Visium CytAssist Gene and Protein Expression Library of Human Glioblastoma, IF, 11mm (FFPE)

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_tissue_image.tif",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_image.tif",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_alignment_file.json",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_molecule_info.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz",
]

##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

subprocess.run("tar -xzf data/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz -C data", shell=True, check=True)
subprocess.run("rm data/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz", shell=True, check=True)
