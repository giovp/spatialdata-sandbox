# G4X Sample Data

Download and convert a sample dataset from the Singular Genomics G4X Platform to SpatialData. The data will include:

- images:
  - `h_and_e`: H&E image in RGB
  - `nuclear`: nuclear stain
  - `eosin`: eosin stain
  - `protein`: multi-channel stack of 20 protein images
- labels:
  - `nuclei`: Nuclei segmentation instance labels
  - `nuclei_exp`: Cell segmentation instance labels
- points:
  - `transcripts`: Single-molecule transcript locations, gene identities, and metadata
- shapes:
  - `nuclei_shapes`: Nuclei segmentation polygons
  - `nuclei_exp`: Cell segmentation polygons
- table:
  - `table`: the cell by gene expression count matrix with cell metadata

### Download
1. Make a folder `g4x_io` and cd into it.
1. Download sample dataset with `download.py`
3. Convert the data into the SpatialData format with `to_zarr.py`
