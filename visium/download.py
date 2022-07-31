#!/usr/bin/env python3

import os
import shutil
from pathlib import Path
import scanpy as sc

from utils import download, unzip

main_dir = Path().resolve() / "data"
print(main_dir)
main_dir.mkdir(parents=True, exist_ok=True)

image_dir = main_dir / "images"
image_dir.mkdir(parents=True, exist_ok=True)

table_dir = main_dir / "tables"
table_dir.mkdir(parents=True, exist_ok=True)


brainfile = main_dir / "mousebrain.zip"
download(
    "https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_visium_wo_cloupe_data.zip",
    brainfile,
    desc="data",
)
unzip(brainfile, main_dir)

for lib in os.scandir(
    main_dir / os.path.join("mouse_brain_visium_wo_cloupe_data", "rawdata")
):
    if not lib.name.startswith("."):
        shutil.copy2(
            main_dir / os.path.join(lib.path, "spatial", "tissue_hires_image.png"),
            main_dir / os.path.join("images", lib.name + ".png"),
        )

        adata = sc.read_visium(main_dir / lib.path)
        adata.write_h5ad(
            main_dir / os.path.join("tables", lib.name + ".h5ad"),
            compression="gzip",
            compression_opts=9,
        )
