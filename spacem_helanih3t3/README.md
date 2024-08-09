### Data

This is a metabolomics dataset from experiments on Hepa and NIH3T3 cells using the [SpaceM](https://doi.org/10.1038/s41592-021-01198-0) method, by [Alexandrov group, EMBL](https://www.embl.org/groups/alexandrov/).

The data consist of the following items:

- coordinate systems:
  Each set of processed images/labels is registered in a corresponding coordinate system with matching prefix.
  (This is because "global" is the default coordinate system of incoming unregistered data and is treated unmutable).
- images:
  - `….pre_maldi`: Microscopy, with `Trans` and `GFP` channels
  - `….post_maldi`: Microscopy after MALDI measurements, with `Trans` and `Dapi` channels
- labels:
  - `….cells`: Segmentation of pre-MALDI images
  - `….ablation_marks`: Segmentation of post-MALDI images
- shapes:
  - `….layout`: Bounding boxes of wells on a slide
  - `….maldi_regions`: Bounding boxes for the MALDI measurements
- tables:
  - `table`:
    - for all annotated elements: `project_id`, `slide_id`, `well_id`, `maldi_region_id`
    - for segmentations:
      - `object_type`, `replicate`, `treatment`
      - scikit-image region properties
      - `X`: MALDI ion intensities

### Download

The dataset is already natively in SpatialData 0.1.2 format.

Download the data with `download.py`, (`write_zarr.py` exists solely for consistency).
