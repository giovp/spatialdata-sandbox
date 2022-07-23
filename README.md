# spatialdata-sandbox

Example datasets for SpatialData development.
## How to add datasets

* Directory per dataset, named after the technology
* Single python executable `download.py` that dump all available data in `data` folder.
* Single python executable `write_zarr.py` that saves all available data in zarr file `data.zarr`.
* Be explicit about where the components are, e.g.
  * Table
  * Points
  * Images
  * Regions
* Have a readme.md with prose