#!/usr/bin/env python3
##
import os
import subprocess
from pathlib import Path

import sys

sys.path.insert(1, os.path.join(sys.path[0], Path(__file__).parent.parent.resolve()))

from utils import download, unzip

##

# from https://www.10xgenomics.com/datasets/multiomic-integration-neuroscience-application-note-visium-for-ffpe-plus-immunofluorescence-alzheimers-disease-mouse-model-brain-coronal-sections-from-one-hemisphere-over-a-time-course-1-standard
# Multiomic Integration Neuroscience Application Note: Visium for FFPE Plus Immunofluorescence Alzheimer's Disease Mouse Model Brain Coronal Sections from One Hemisphere Over a Time Course

urls = [
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.1/VisiumFFPE_Mouse_Brain_Alzheimers_AppNote/VisiumFFPE_Mouse_Brain_Alzheimers_AppNote_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.1/VisiumFFPE_Mouse_Brain_Alzheimers_AppNote/VisiumFFPE_Mouse_Brain_Alzheimers_AppNote_spatial.tar.gz",
]

##
os.makedirs("data", exist_ok=True)
for url in urls:
    print(url)
    name = Path(url).name
    download(url, os.path.join("data", name), name)

subprocess.run("tar -xzf data/VisiumFFPE_Mouse_Brain_Alzheimers_AppNote_spatial.tar.gz -C data", shell=True, check=True)
subprocess.run("rm data/VisiumFFPE_Mouse_Brain_Alzheimers_AppNote_spatial.tar.gz", shell=True, check=True)
