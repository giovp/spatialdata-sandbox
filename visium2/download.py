#!/usr/bin/env python3
##
import os
import shutil
from pathlib import Path
import scanpy as sc

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

##

# from https://www.10xgenomics.com/resources/datasets/adult-mouse-olfactory-bulb-1-standard-1
# Adult Mouse Olfactory Bulb
# Spatial Gene Expression Dataset by Space Ranger 2.0.0

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_molecule_info.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz",
]

##
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

os.system("tar -xzf data/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz -C data")
# scanpy is expecting an old version of the spaceranger output, so we termporarily rename the files until we refactor
# the parser
os.system("rm data/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz")
os.system(
    "mv data/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5 data/filtered_feature_bc_matrix.h5"
)
os.system("mv data/Visium_Mouse_Olfactory_Bulb_molecule_info.h5 data/molecule_info.h5")
os.system('mv data/spatial/tissue_positions.csv data/spatial/tissue_positions_list.csv')
