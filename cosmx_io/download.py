##
import os
from pathlib import Path
import subprocess


nanostring_dir = Path().resolve() / "data"
nanostring_dir.mkdir(parents=True, exist_ok=True)
assert nanostring_dir.exists()
path = nanostring_dir / "data_lung5_rep2"

if not path.exists():
    url = "https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep2/Lung5_Rep2+SMI+Flat+data.tar.gz"
    ugly_name = "Lung5_Rep2/Lung5_Rep2-Flat_files_and_images/"
    subprocess.run(f"wget -P {nanostring_dir} {url}", shell=True, check=True)
    subprocess.run(
        f"tar -xzf {nanostring_dir/'Lung5_Rep2+SMI+Flat+data.tar.gz'} -C {nanostring_dir}",
        shell=True,
        check=True,
    )
    subprocess.run(f"mv {nanostring_dir/ugly_name} {path}", shell=True, check=True)
    subprocess.run(f"rm -r {nanostring_dir/'Lung5_Rep2'}", shell=True, check=True)
    subprocess.run(
        f"rm {nanostring_dir/'Lung5_Rep2+SMI+Flat+data.tar.gz'}", shell=True, check=True
    )

counts_file = "Lung5_Rep2_exprMat_file.csv"
meta_file = "Lung5_Rep2_metadata_file.csv"
fov_file = "Lung5_Rep2_fov_positions_file.csv"
transcripts_file = "Lung5_Rep2_tx_file.csv"
