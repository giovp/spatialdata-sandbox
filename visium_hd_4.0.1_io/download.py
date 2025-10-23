#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download

##

# from https://www.10xgenomics.com/datasets/visium-hd-three-prime-mouse-brain-fresh-frozen
# Visium HD 3' Gene Expression Library, Mouse Brain (Fresh Frozen)

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_barcode_mappings.parquet",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_binned_outputs.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_cloupe_008um.cloupe",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_cloupe_cell.cloupe",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_feature_slice.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_metrics_summary.csv",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_molecule_info.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_segmented_outputs.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_spatial.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_web_summary.html",
]


##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

print("Unpacking files...")
# Unpack tar.gz files
for fname in os.listdir("data"):
    if fname.endswith(".tar.gz"):
        tar_path = os.path.join("data", fname)
        subprocess.run(["tar", "-xzf", tar_path, "-C", "data"], check=True)

print("Download and unpack complete.")
