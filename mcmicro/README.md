# MCMICRO dataset (CODEX+cycif+mIHC)

## Dataset overview
This dataset is a multimodal dataset consisting out of whole slide images
of sequential sections of a human tonsil specimen 
([MCMICRO](https://www.nature.com/articles/s41592-021-01308-y)). Only the CODEX dataset is downloaded. Each sequential
section is imaged using a different imaging technology. Raw unstitched data is not available. If useful it can be requested from Denis Shapiro. For an overview of the 
sections in the dataset, see image (b-d) below. The size of the entire dataset is 116.5GB. 
<img src='https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41592-021-01308-y/MediaObjects/41592_2021_1308_Fig2_HTML.png?as=webp' alt='MCMICRO' title='MCMICRO sections'/>

## Requirements
Create an account at [Synapse](https://www.synapse.org/#!RegisterAccount:0) and run the following command

`pip install synapseclient squidpy`

The synapseclient is required to download the data.

## Downloading the data
The data can be downloaded using the synapse cli. However, this is not adviced. Instead use *synapse_cli.py* as follows:
 
`python synapse_cli.py -u <USERNAME> -p <PASSWORD>`

This will save the datasets in the data directory as follows:

```bash
.
├── data
│   └── WSI_tonsil
│       └── CODEX
│           ├── CODEX.h5ad
│           ├── markers.csv
│           ├── probability-maps
│           │   ├── SYNAPSE_METADATA_MANIFEST.tsv
│           │   └── unmicst
│           │       ├── PilotTonsil_5_z08_Probabilities_0.tif
│           │       └── SYNAPSE_METADATA_MANIFEST.tsv
│           ├── qc
│           │   ├── params.yml
│           │   ├── provenance
│           │   │   ├── quantification.log
│           │   │   ├── quantification.sh
│           │   │   ├── s3seg.log
│           │   │   ├── s3seg.sh
│           │   │   ├── SYNAPSE_METADATA_MANIFEST.tsv
│           │   │   ├── unmicst.log
│           │   │   └── unmicst.sh
│           │   └── SYNAPSE_METADATA_MANIFEST.tsv
│           ├── quantification
│           │   ├── SYNAPSE_METADATA_MANIFEST.tsv
│           │   └── unmicst-PilotTonsil_5_z08.csv
│           ├── registration
│           │   ├── PilotTonsil_5_z08.ome.tif
│           │   └── SYNAPSE_METADATA_MANIFEST.tsv
│           ├── segmentation
│           │   ├── SYNAPSE_METADATA_MANIFEST.tsv
│           │   └── unmicst-PilotTonsil_5_z08
│           │       ├── cellMask.tif
│           │       ├── cellRingMask.tif
│           │       ├── cytoMask.tif
│           │       ├── cytoRingMask.tif
│           │       ├── nucleiMask.tif
│           │       ├── nucleiRingMask.tif
│           │       └── SYNAPSE_METADATA_MANIFEST.tsv
│           └── SYNAPSE_METADATA_MANIFEST.tsv
```

Afterwards run *mcmicro.py* to store the h5ad file in each TNP_pilot directory.
