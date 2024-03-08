#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

##

# from https://www.10xgenomics.com/datasets/adult-human-brain-1-cerebral-cortex-unknown-orientation-stains-anti-gfap-anti-nfh-1-standard-1-1-0
# Adult Human Brain 1 (Cerebral Cortex, Unknown orientation).

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_1/V1_Human_Brain_Section_1_image.tif",
    # "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_1/V1_Human_Brain_Section_1_alignment_file.json",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_1/V1_Human_Brain_Section_1_molecule_info.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_1/V1_Human_Brain_Section_1_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_1/V1_Human_Brain_Section_1_spatial.tar.gz",
]

##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

subprocess.run("tar -xzf data/V1_Human_Brain_Section_1_spatial.tar.gz -C data", shell=True, check=True)
subprocess.run("rm data/V1_Human_Brain_Section_1_spatial.tar.gz", shell=True, check=True)
