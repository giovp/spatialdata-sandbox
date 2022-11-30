### Requirements
`pip install anndata datashader imageio`

### Data
MERFISH M. Musculus (VISp), Allen Institute prototype MERFISH pipeline. Downloaded from [spacetx-spacejam] (https://github.com/spacetx-spacejam/data).

Raw data saved in `raw/`

Processed data saved in
- `processed/single_molecule.h5ad` (AnnData with single molecule points coordinates and annotation)
- `processed/cells.h5ad` (AnnData with single cell regions and gene expression)
- `processed/anatomical.geojson` (cortex layer anatomical regions annotation)
- `processed/image.png` (rasterized version of the point data)
- `processed/image_transorm.json` translation and scaling parameters for aligning the image with the rest of the coordinates

Set in the code `PLOT = True` for some basic visualization of the data.

Download the data with `download.py`, then save it to the OME-NGFF format with `write_zarr.py`.
