
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import spatialdata_plot

    path_read = "/Users/macbook/ssd/biodata/seqfish/instrument 2 official"
    sdata = seqfish(
        path=path_read,
        rois=[1],
        cells_as_circles=True,
    )

    ##
    axes = plt.subplots(1, 3, figsize=(15, 5))[1]
    gene_name = "Arg1"

    # plot gene expression over segmentation labels
    (
        sdata.pl.render_images("Roi1_DAPI", cmap="gray")
        .pl.render_labels("Roi1_Segmentation", color=gene_name)
        .pl.render_points("Roi1_TranscriptList", color="name", groups=gene_name, palette="orange")
        .pl.show(
            title="Plotting cells as labels", coordinate_systems="global", figsize=(10, 5), ax=axes[0]
        )
    )

    # plot gene expression over boundaries shapes
    sdata["table_Roi1"].obs["region"] = "Roi1_Boundaries"
    sdata.set_table_annotates_spatialelement(
        table_name="table_Roi1", region="Roi1_Boundaries", region_key="region", instance_key="instance_id"
    )

    (
        sdata.pl.render_images("Roi1_DAPI", cmap="gray")
        .pl.render_shapes("Roi1_Boundaries", color=gene_name)
        .pl.render_points("Roi1_TranscriptList", color="name", groups=gene_name, palette="orange")
        .pl.show(
            title="Plotting cells as shapes", coordinate_systems="global", figsize=(10, 5), ax=axes[1]
        )
    )

    # plot gene expression over cell circles
    sdata["table_Roi1"].obs["region"] = "Roi1_CellCoordinates"
    sdata.set_table_annotates_spatialelement(
        table_name="table_Roi1", region="Roi1_CellCoordinates", region_key="region", instance_key="instance_id"
    )
    (
        sdata.pl.render_images("Roi1_DAPI", cmap="gray")
        .pl.render_shapes("Roi1_CellCoordinates", color=gene_name)
        .pl.render_points("Roi1_TranscriptList", color="name", groups=gene_name, palette="orange")
        .pl.show(
            title="Plotting cells as circles", coordinate_systems="global", figsize=(10, 5), ax=axes[2]
        )
    )

    plt.suptitle(f'{gene_name} expression over DAPI image')
    plt.tight_layout()
    plt.show()
    ##

    from napari_spatialdata import Interactive

    Interactive(sdata)

