##
from pathlib import Path
from spatialdata_io import stereoseq
import shutil
import spatialdata as sd

##
path = Path().resolve()
if not str(path).endswith("stereoseq_io"):
    path /= "stereoseq_io"
    assert path.exists()
path_write = path / "data.zarr"

WRITE = True
# WRITE = False
if WRITE:
    ##
    # you can use symlinks to make the data available
    # C://STT1_image_alignment//pipeline output//result"
    dataset = Path("data/STT1_image_alignment/pipeline output/result")

    path_read = path / dataset
    sdata = stereoseq(path_read)

    ##
    print("writing the data... ", end="")
    if path_write.exists():
        shutil.rmtree(path_write)
    sdata.write(path_write)
    print("done")

##
print(f'view with "python -m napari_spatialdata view data.zarr"')
sdata = sd.SpatialData.read(path_write)
print("read")
print(sdata)

# ##
# adata = sdata['bin10_table']
# import scanpy as sc
# adata.X = adata.X.astype(float)
# adata.X = sc.pp.log1p(adata.X)
#
# sc.tl.pca(adata, svd_solver="arpack")
# sc.pl.pca_variance_ratio(adata)
#
# ##
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
#
# ##
# sc.tl.leiden(
#     adata,
#     resolution=0.9,
#     random_state=0,
#     flavor="igraph",
#     n_iterations=2,
#     directed=False,
# )
# ##
# sdata['bin10_table_processed'] = adata

##
from napari_spatialdata.constants import config
config.POINT_THRESHOLD = 1000000000
from napari_spatialdata import Interactive
Interactive(sdata)
