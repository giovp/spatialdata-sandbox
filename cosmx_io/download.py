##
import os
import numpy as np
import pandas as pd
import scanpy as sc
import time

from anndata import AnnData
from pathlib import Path
from scipy.sparse import csr_matrix

nanostring_dir = Path().resolve() / "data"
nanostring_dir.mkdir(parents=True, exist_ok=True)
assert nanostring_dir.exists()
path = nanostring_dir / "data_lung5_rep2"
# if pyarrow is available, better to fail now, than after minutes of downloads
import pyarrow

# to prevent the unused import to be accidentally removed by PyCharm's "optimize imports"
_ = pyarrow.__path__

if not path.exists():
    url = "https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep2/Lung5_Rep2+SMI+Flat+data.tar.gz"
    ugly_name = "Lung5_Rep2/Lung5_Rep2-Flat_files_and_images/"
    os.system(f"wget -P {nanostring_dir} {url}")
    os.system(
        f"tar -xzf {nanostring_dir/'Lung5_Rep2+SMI+Flat+data.tar.gz'} -C {nanostring_dir}"
    )
    os.system(f"mv {nanostring_dir/ugly_name} {path}")
    os.system(f"rm -r {nanostring_dir/'Lung5_Rep2'}")
    os.system(f"rm {nanostring_dir/'Lung5_Rep2+SMI+Flat+data.tar.gz'}")

counts_file = "Lung5_Rep2_exprMat_file.csv"
meta_file = "Lung5_Rep2_metadata_file.csv"
fov_file = "Lung5_Rep2_fov_positions_file.csv"
transcripts_file = "Lung5_Rep2_tx_file.csv"
