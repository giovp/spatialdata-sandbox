##
from napari_spatialdata import Interactive
import spatialdata as sd

sdata = sd.read_zarr("data.zarr")
print(sdata)
##
import pandas as pd
from pathlib import Path

parent_dir = Path(__file__).parent
bins_size = "002"

# add clusters
# f = parent_dir / f'data/square_{bins_size}um/analysis/clustering/gene_expression_graphclust/clusters.csv'
# assert f.exists()

# df = pd.read_csv(f)
# df
# df["Cluster"] = df["Cluster"].astype("category")
# df.set_index("Barcode", inplace=True)
# sdata[f"square_{bins_size}um"].obs["Cluster"] = df["Cluster"]

##

sdata[f'square_{bins_size}um'].X = sdata[f'square_{bins_size}um'].X.tocsc()
bins = sd.rasterize_bins(
    sdata,
    bins=f"Visium_HD_Human_Lung_Cancer_Fixed_Frozen_square_{bins_size}um",
    table_name=f'square_{bins_size}um',
    col_key="array_col",
    row_key="array_row",
)
##

sdata[f'rasterized_{bins_size}um'] = bins
##
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from napari.utils.colormaps import Colormap

import matplotlib.pyplot as plt

cmap = plt.get_cmap("viridis")
cmap_zero = set_zero_in_cmap_to_transparent(cmap)
napari_cmap = Colormap(cmap_zero.colors, 'viridis_zero_transparent')

interactive = Interactive(sdata, headless=True)
interactive.add_element(f'rasterized_{bins_size}um', 'global', view_element_system=True)
interactive.get_layer(f'rasterized_{bins_size}um').colormap = napari_cmap
interactive.run()