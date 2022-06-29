#!/usr/bin/env python3

import os
import shutil

import scanpy as sc

from utils import download, unzip

os.makedirs("images")

os.makedirs("tables")
brainfile = "mousebrain.zip"
download(
    "https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_visium_wo_cloupe_data.zip",
    brainfile,
    desc="data",
)
unzip(brainfile)

for lib in os.scandir(os.path.join("mouse_brain_visium_wo_cloupe_data", "rawdata")):
    if not lib.name.startswith("."):
        shutil.copy2(os.path.join(lib.path, "spatial", "tissue_hires_image.png"), os.path.join("images", lib.name + ".png"))

        adata = sc.read_visium(lib.path)
        adata.write_h5ad(os.path.join("tables", lib.name + ".h5ad"), compression="gzip", compression_opts=9)
