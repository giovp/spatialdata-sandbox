#!/usr/bin/env python3

import os
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

main_dir = Path().resolve() / "data"
print(main_dir)
main_dir.mkdir(parents=True, exist_ok=True)


brainfile = main_dir / "mousebrain.zip"
download(
    "https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_visium_wo_cloupe_data.zip",
    brainfile,
    desc="data",
)
unzip(brainfile, main_dir)