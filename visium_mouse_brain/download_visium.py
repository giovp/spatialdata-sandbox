#!/usr/bin/env python3

import os
import urllib.request
import zipfile
import shutil

from tqdm import tqdm
import scanpy as sc


class TqdmDownload(tqdm):
    def __init__(self, *args, **kwargs):
        kwargs = dict(kwargs)
        kwargs.update({"unit": "B", "unit_scale": True, "unit_divisor": 1024})
        super().__init__(*args, **kwargs)

    def update_to(self, nblocks=1, blocksize=1, total=-1):
        self.total = total
        self.update(nblocks * blocksize - self.n)


def download(url, outfile, desc):
    with TqdmDownload(desc="downloading " + desc) as t:
        urllib.request.urlretrieve(url, outfile, t.update_to)


def unzip(file, outdir=None, files=None, rm=True):
    if outdir is None:
        outdir = os.getcwd()
    else:
        os.makedirs(outdir, exist_ok=True)
    zfile = zipfile.ZipFile(file)
    if files is not None:
        for f in files:
            zfile.extract(f, outdir)
    else:
        zfile.extractall(outdir)
    zfile.close()
    if rm:
        os.unlink(file)


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
