# run the following code in a cell at the end of the densenet.ipynb documentation notebook


x = np.array([13694.0, 13889.0, 13889.0, 13694.0, 13694.0])
y = np.array([13984.0, 13984.0, 14162.0, 14162.0, 13984.0])

small_sdata = merged.query.bounding_box(
    axes=("x", "y"),
    min_coordinate=[np.min(x), np.min(y)],
    max_coordinate=[np.max(x), np.max(y)],
    target_coordinate_system="aligned",
)
small_sdata

Interactive(small_sdata)

small_dataset = ImageTilesDataset(
    sdata=small_sdata,
    regions_to_images={"cell_boundaries": "CytAssist_FFPE_Human_Breast_Cancer_full_image"},
    tile_dim_in_units=3 * xenium_circles_diameter,
    tile_dim_in_pixels=32,
    target_coordinate_system="aligned",
    transform=None,
)

small_dataset[0]

import matplotlib.pyplot as plt
import spatialdata as sd
import spatialdata_plot

n = len(small_dataset)
axes = plt.subplots(1, n, figsize=(15, 3))[1]
for sdata_tile, i in zip(small_dataset, range(n)):
    # BUG: we need to explicitly remove the coordinate system global if we want to combine images and shapes plots into a single subplot
    if "global" in get_transformation(sdata_tile["cell_boundaries"], get_all=True):
        sd.transformations.remove_transformation(sdata_tile["cell_boundaries"], "global")
    sdata_tile.pl.render_images().pl.render_shapes(outline_width=3.0, outline=True, fill_alpha=0.0).pl.show(ax=axes[i])
