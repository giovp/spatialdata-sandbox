##
import re
import os
import numpy as np
import spatialdata as sd
from napari_spatialdata import Interactive
from dask_image.imread import imread
import pandas as pd
from anndata import AnnData
import anndata as ad

f = "data/Kidney Data 4 Download"

##
CELL_COORDINATES = "CellCoordinates"
SECTION = "Section"
CXG = "CxG"
DAPI = "DAPI"


def find_library_id_and_sections(data_root):
    r_cell_coordinates = rf"^([A-Za-z0-9_]*?)_{CELL_COORDINATES}_{SECTION}([0-9]+).csv$"
    files = os.listdir(data_root)
    patterns = set()
    sections = set()
    for file in files:
        m = re.match(r_cell_coordinates, file)
        if m:
            patterns.add(m.group(1))
            sections.add(int(m.group(2)))
    assert len(patterns) == 1
    return patterns.pop(), sections


def seqfish(data_root):
    library_id, sections = find_library_id_and_sections(data_root)
    # sections.remove(2)
    # sections.remove(3)
    print(library_id, sections)
    images = {}
    shapes = {}
    tables = {}

    for section in sections:
        cell_coordinates_file = (
            f"{library_id}_{CELL_COORDINATES}_{SECTION}{section}.csv"
        )
        cxg_file = f"{library_id}_{CXG}_{SECTION}{section}.csv"
        image_file = f"{library_id}_{DAPI}_{SECTION}{section}.tiff"
        assert os.path.isfile(os.path.join(data_root, cell_coordinates_file))
        assert os.path.isfile(os.path.join(data_root, cxg_file))
        assert os.path.isfile(os.path.join(data_root, image_file))
        coords = pd.read_csv(os.path.join(data_root, cell_coordinates_file))
        xy = coords[["center_x", "center_y"]].values
        radius = np.sqrt(coords["area"] / np.pi)
        label = coords["label"].values
        cxg = pd.read_csv(os.path.join(data_root, cxg_file))
        im = imread(os.path.join(data_root, image_file))
        images[f"image_Section{section}"] = sd.models.Image2DModel.parse(
            im,
            dims=["c", "y", "x"],
            scale_factors=[2, 2, 2, 2],
            transformations={f"Section{section}": sd.transformations.Identity()},
        )
        shapes[f"cells_Section{section}"] = sd.models.ShapesModel.parse(
            xy,
            geometry=0,
            radius=radius,
            transformations={f"Section{section}": sd.transformations.Identity()},
            index=label,
        )
        table = AnnData(
            cxg.iloc[:, 1:],
            obs=pd.DataFrame(
                {"region": f"cells_Section{section}", "instance_id": cxg.iloc[:, 0]}
            ),
        )
        tables[f"Section{section}"] = table

    merged = ad.concat(tables)
    table = sd.models.TableModel.parse(
        merged,
        region=[f"cells_Section{section}" for section in sections],
        region_key="region",
        instance_key="instance_id",
    )

    return sd.SpatialData(images=images, shapes=shapes, table=table)


if __name__ == "__main__":
    sdata = seqfish(f)
    print(sdata)
    sdata.write('data.zarr')
    # interactive = Interactive(sdata)
    # interactive.run()

