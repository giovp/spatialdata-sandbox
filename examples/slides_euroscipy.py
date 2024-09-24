##
import spatialdata as sd
import matplotlib.pyplot as plt
from napari_spatialdata import Interactive
import spatialdata_plot

xenium_path = "xenium_rep1_io/data_aligned.zarr"

sdata = sd.read_zarr(xenium_path)

##

min_coordinate = [11750, 32549]
max_coordinate = [13880, 35114]

sdata_small = sdata.query.bounding_box(
    axes=("y", "x"),
    min_coordinate=min_coordinate,
    max_coordinate=max_coordinate,
    target_coordinate_system="aligned",
)
sdata_small["table"].obs["region"] = "cell_boundaries"
sdata_small.set_table_annotates_spatialelement(
    table_name="table",
    region="cell_boundaries",
    region_key="region",
    instance_key="cell_id",
)

# from napari_spatialdata import Interactive
##

sdata_small.write("sdata_small.zarr", overwrite=True)
sdata_small = sd.read_zarr("sdata_small.zarr")

##
# bug: cannot set cmap to gray for multiscale image
sdata_small["morphology_focus_scale0"] = sdata_small["morphology_focus"]["scale0"][
    "image"
]

##
# sdata_small.pl.render_images('morphology_focus_scale0', cmap='gray').pl.render_points('transcripts', method='datashader').pl.show(coordinate_systems='aligned')
# plt.show()
transcripts = sdata_small["transcripts"].compute()
categories = transcripts["feature_name"].cat.categories[:20]
some_points = transcripts[transcripts["feature_name"].isin(categories)]
transformations = sd.transformations.get_transformation(
    sdata_small["transcripts"], get_all=True
)
some_points["feature_name"] = some_points["feature_name"].cat.remove_unused_categories()
parsed = sd.models.PointsModel.parse(some_points, transformations=transformations)
sdata_small["some_points"] = parsed

##

# bug: with spatialdata-plot when instances and table rows don't match
filtered_shapes, filtered_table = sd.join_spatialelement_table(
    sdata_small,
    spatial_element_names="cell_boundaries",
    table_name="table",
    how="inner",
)
sdata_small["filtered_shapes"] = filtered_shapes["cell_boundaries"]
sdata_small["filtered_table"] = filtered_table
filtered_table.obs["region"] = "filtered_shapes"
filtered_table.uns["spatialdata_attrs"]["region"] = "filtered_shapes"

##
# bug: outline doesn't color by celltype_major
plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata_small.pl.render_images("morphology_focus_scale0", cmap="gray").pl.render_points(
    "some_points", color="feature_name", method="matplotlib"
).pl.render_shapes(
    "filtered_shapes",
    color="celltype_major",
    fill_alpha=0.0,
    # outline=True,
    table_name="filtered_table",
).pl.show(
    coordinate_systems="aligned", ax=ax
)
plt.show()

# Interactive(sdata_small)
