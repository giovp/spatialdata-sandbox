
import json
import os
import shutil
from pathlib import Path
import sys
sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download

# downloads a codex minimal example
print(f"os.getcwd() = {os.getcwd()}")
data_dir = Path().resolve() / "data"
data_dir.mkdir(parents=True, exist_ok=True)
assert data_dir.exists()
raw = data_dir / "raw"
raw.mkdir(parents=True, exist_ok=True)
processed = data_dir / "processed"
processed.mkdir(parents=True, exist_ok=True)



BASE_URL = 'https://www.huber.embl.de/users/harald/codex/'
FILES = ['annotations.csv' , 'expression.csv',  'markerlist.csv', 'lymphnode.tif', 'segmentation.tif',]

def my_download(url):
    # aria2c asks for relative paths
    download(
        url,
        (raw / os.path.basename(url)).relative_to(Path().resolve()),
        desc="data",
    )


for f in FILES:
    if os.path.exists(raw / f):
        print(f"Skipping {f}")
        continue
    
    my_download(f'{BASE_URL}{f}')
