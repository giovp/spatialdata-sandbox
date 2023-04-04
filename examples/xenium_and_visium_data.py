##
from anndata import AnnData
import torch.nn.functional as F
import numpy as np
from spatial_image import SpatialImage
import spatialdata as sd
from napari_spatialdata import Interactive
from spatialdata import SpatialData
from spatialdata.transformations import get_transformation
from spatialdata.dataloader.datasets import ImageTilesDataset
from tqdm import tqdm
from spatialdata import transform
from typing import Dict
import os
import scanpy as sc

from monai.networks.nets import DenseNet121
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor
from pytorch_lightning.callbacks.progress import TQDMProgressBar
from pytorch_lightning.loggers import CSVLogger
from torch.nn import CrossEntropyLoss
from torch.optim import Adam
import torchvision
import torch
from torch.utils.data import DataLoader
from pytorch_lightning import LightningDataModule

##
SMALL_DATASET = False


def load_data() -> tuple[SpatialData, SpatialData, SpatialData]:
    # loading the Xenium and Visium SpatialData objects
    print("current working directory:", os.getcwd())
    SPATIALDATA_SANDBOX_PATH = "."
    GENERATED_DATA_PATH = os.path.join(
        SPATIALDATA_SANDBOX_PATH, "generated_data/xenium_visium_integration"
    )
    assert os.path.isdir(
        GENERATED_DATA_PATH
    ), f"{GENERATED_DATA_PATH} not found, please use symlinks to make it available"

    XENIUM_SDATA_PATH = os.path.join(
        SPATIALDATA_SANDBOX_PATH, "xenium_rep1_io/data.zarr"
    )
    VISIUM_SDATA_PATH = os.path.join(
        SPATIALDATA_SANDBOX_PATH, "visium_associated_xenium_io/data.zarr"
    )

    assert os.path.isdir(XENIUM_SDATA_PATH)
    assert os.path.isdir(VISIUM_SDATA_PATH)

    xenium_sdata = sd.read_zarr(XENIUM_SDATA_PATH)
    visium_sdata = sd.read_zarr(VISIUM_SDATA_PATH)

    VISIUM_CLONAL_DATA = os.path.join(GENERATED_DATA_PATH, "visium_copyKat.h5ad")
    assert os.path.isfile(VISIUM_CLONAL_DATA)
    adata = sc.read_h5ad(VISIUM_CLONAL_DATA)

    ##
    # create a dataset with tiles, one for each Visium circle, but around the Xenium image
    merged = sd.SpatialData(
        images={
            "xenium": xenium_sdata.images["morphology_mip"],
            "visium": visium_sdata.images[
                "CytAssist_FFPE_Human_Breast_Cancer_full_image"
            ],
        },
        shapes={
            "CytAssist_FFPE_Human_Breast_Cancer": visium_sdata.shapes[
                "CytAssist_FFPE_Human_Breast_Cancer"
            ]
        },
        table=visium_sdata.table,
    )

    s = adata.obs["clone"]
    s.index = merged.table.obs.index
    visium_sdata.table.obs["clone"] = s

    if SMALL_DATASET:
        # subsetting the data for debug/development purposes
        min_coordinate = [12790, 12194]
        max_coordinate = [15100, 14221]
        merged = merged.query.bounding_box(
            min_coordinate=min_coordinate,
            max_coordinate=max_coordinate,
            axes=["y", "x"],
            target_coordinate_system="aligned",
        )
    return merged, visium_sdata, xenium_sdata


def _get_clonal_classes(adata: AnnData) -> list[str]:
    categories = adata.obs["clone"].cat.categories.tolist()
    return categories


class ClonalTileDataset(torch.utils.data.Dataset):
    def __init__(self, tile_dataset: ImageTilesDataset):
        self.tile_dataset = tile_dataset
        self.categories = _get_clonal_classes(self.tile_dataset.sdata.table)

    def __len__(self):
        return len(self.tile_dataset)

    def __getitem__(self, idx):
        tile, row = self.tile_dataset[idx]
        clonal_information = row.obs["clone"].values[0]
        clonal_information = self.categories.index(clonal_information)
        clonal_one_hot = F.one_hot(
            torch.tensor(clonal_information), num_classes=len(self.categories)
        )
        tile_numpy = tile.data.compute().astype(np.float32)
        return tile_numpy, clonal_one_hot.float()


def get_dataset() -> ImageTilesDataset:
    merged, xenium_sdata, visium_sdata = load_data()

    circles = merged["CytAssist_FFPE_Human_Breast_Cancer"]
    t = get_transformation(circles, "aligned")

    transformed_circles = transform(circles, t)
    visium_circle_diameter = 2 * transformed_circles.radius.iloc[0]
    dataset = ImageTilesDataset(
        sdata=merged,
        regions_to_images={"CytAssist_FFPE_Human_Breast_Cancer": "xenium"},
        tile_dim_in_units=visium_circle_diameter,
        tile_dim_in_pixels=32,
        target_coordinate_system="aligned",
    )
    ##
    # visualization (only if using the small dataset, otherwise the code will be too slow)
    if SMALL_DATASET:
        tiles = []
        for i, tile in enumerate(tqdm(dataset)):
            tiles.append(tile)

        tiles_sdata = SpatialData(
            images={
                "full": visium_sdata.images[
                    "CytAssist_FFPE_Human_Breast_Cancer_full_image"
                ]
            }
            | {
                f"tile_{region_name}_{region_index}": tile
                for (tile, region_name, region_index) in tiles
            }
        )

        Interactive([tiles_sdata, merged])
    ##

    clonal_tile_dataset = ClonalTileDataset(dataset)
    return clonal_tile_dataset


##


class TilesDataModule(LightningDataModule):
    def __init__(
        self,
        batch_size: int,
        num_workers: int,
        train_transform: torchvision.transforms.Compose,
        test_transform: torchvision.transforms.Compose,
    ):
        super().__init__()
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.train_transform = train_transform
        self.test_transform = test_transform
        self.train = None
        self.val = None
        self.test = None

    def setup(self, stage=None):
        self.dataset = get_dataset()
        n_train = int(len(self.dataset) * 0.7)
        n_val = int(len(self.dataset) * 0.2)
        n_test = len(self.dataset) - n_train - n_val
        self.train, self.val, self.test = torch.utils.data.random_split(
            self.dataset,
            [n_train, n_val, n_test],
            generator=torch.Generator().manual_seed(42),
        )
        self.train.dataset.tile_dataset.transform = self.train_transform
        self.val.dataset.tile_dataset.transform = self.test_transform
        self.test.dataset.tile_dataset.transform = self.test_transform

    def train_dataloader(self):
        return DataLoader(
            self.train,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            shuffle=True,
        )

    def val_dataloader(self):
        return DataLoader(
            self.val,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            shuffle=False,
        )

    def test_dataloader(self):
        return DataLoader(
            self.test,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            shuffle=False,
        )

    def predict_dataloader(self):
        return DataLoader(
            self.dataset,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            shuffle=False,
        )


class DenseNetModel(pl.LightningModule):
    def __init__(self, learning_rate: float, in_channels: int, num_classes: int):
        super().__init__()

        # store hyperparameters
        self.save_hyperparameters()

        self.loss_function = CrossEntropyLoss()

        # make the model
        self.model = DenseNet121(
            spatial_dims=2, in_channels=in_channels, out_channels=num_classes
        )

    def forward(self, x) -> torch.Tensor:
        return self.model(x)

    def _compute_loss_from_batch(
        self, batch: Dict[str, torch.Tensor], batch_idx: int
    ) -> float:
        inputs = batch[0]
        labels = batch[1]

        outputs = self.model(inputs)
        return self.loss_function(outputs, labels)

    def training_step(
        self, batch: Dict[str, torch.Tensor], batch_idx: int
    ) -> Dict[str, float]:
        # compute the loss
        loss = self._compute_loss_from_batch(batch=batch, batch_idx=batch_idx)

        # perform logging
        self.log("training_loss", loss, batch_size=len(batch[0]))

        return {"loss": loss}

    def validation_step(self, batch: Dict[str, torch.Tensor], batch_idx: int) -> float:
        loss = self._compute_loss_from_batch(batch=batch, batch_idx=batch_idx)

        imgs, labels = batch
        acc = self.compute_accuracy(imgs, labels)
        # By default logs it per epoch (weighted average over batches), and returns it afterwards
        self.log("test_acc", acc)

        return loss

    def test_step(self, batch, batch_idx):
        imgs, labels = batch
        acc = self.compute_accuracy(imgs, labels)
        # By default logs it per epoch (weighted average over batches), and returns it afterwards
        self.log("test_acc", acc)

    def predict_step(self, batch, batch_idx: int, dataloader_idx: int = 0):
        imgs, labels = batch
        preds = self.model(imgs).argmax(dim=-1)
        return preds

    def compute_accuracy(self, imgs, labels):
        preds = self.model(imgs).argmax(dim=-1)
        labels_value = torch.argmax(labels, dim=-1)
        acc = (labels_value == preds).float().mean()
        return acc

    def configure_optimizers(self) -> Adam:
        return Adam(self.model.parameters(), lr=self.hparams.learning_rate)
