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
import torch.multiprocessing as mp


##
# to fix a mysterious bug with dataloaders in the gpu machine
mp.set_start_method("spawn", force=True)  # You can also try 'fork' or 'forkserver'

# SMALL_DATASET = False
SMALL_DATASET = True


def load_raw_data() -> tuple[SpatialData, SpatialData]:
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
    s = adata.obs["clone"]
    s.index = visium_sdata.table.obs.index
    visium_sdata.table.obs["clone"] = s

    return xenium_sdata, visium_sdata


class TileDataset(torch.utils.data.Dataset):
    CATEGORICAL_COLUMN = "celltype_major"

    def __init__(self, transform = None):
        self.transform = transform
        self.dataset0, self.dataset1 = self._init_tiles()
        self.categories = self.dataset0.sdata.table.obs[
            self.CATEGORICAL_COLUMN
        ].cat.categories.tolist()

    def _init_tiles(self) -> ImageTilesDataset:
        xenium_sdata, visium_sdata = load_raw_data()

        # create a dataset with tiles, one for each Xenium cell, but around the Visium H&E image
        merged = sd.SpatialData(
            images={
                "morphology_mip": xenium_sdata.images["morphology_mip"],
                "CytAssist_FFPE_Human_Breast_Cancer_full_image": visium_sdata.images[
                    "CytAssist_FFPE_Human_Breast_Cancer_full_image"
                ],
            },
            shapes={"cell_circles": xenium_sdata.shapes["cell_circles"]},
            table=xenium_sdata.table,
        )
        if SMALL_DATASET:
            merged = self._crop_dataset(merged)

        circles = merged["cell_circles"]
        t = get_transformation(circles, "aligned")

        transformed_circles = transform(circles, t)
        xenium_circles_diameter = 2 * np.mean(transformed_circles.radius)

        dataset0 = ImageTilesDataset(
            sdata=merged,
            regions_to_images={
                "cell_circles": "CytAssist_FFPE_Human_Breast_Cancer_full_image"
            },
            tile_dim_in_units=4 * xenium_circles_diameter,
            tile_dim_in_pixels=32,
            target_coordinate_system="aligned",
        )

        dataset1 = ImageTilesDataset(
            sdata=merged,
            regions_to_images={"cell_circles": "morphology_mip"},
            tile_dim_in_units=2 * xenium_circles_diameter,
            tile_dim_in_pixels=32,
            target_coordinate_system="aligned",
        )

        if SMALL_DATASET and False:
            self._napari_visualization(merged, dataset0)
            self._napari_visualization(merged, dataset1)
            os._exit(0)

        visium_image = merged['CytAssist_FFPE_Human_Breast_Cancer_full_image']['scale0']['CytAssist_FFPE_Human_Breast_Cancer_full_image'].data.compute()
        xenium_image = merged['morphology_mip']['scale0']['morphology_mip'].data.compute()
        print('computing standard deviation and mean...', end='')
        self.visium_mean = np.mean(visium_image.reshape(visium_image.shape[0], -1), axis=-1)
        self.visium_std = np.std(visium_image.reshape(visium_image.shape[0], -1), axis=-1)
        self.xenium_mean = np.mean(xenium_image)
        self.xenium_std = np.std(xenium_image)
        print('done')
        return dataset0, dataset1

    @staticmethod
    def _crop_dataset(merged: SpatialData) -> SpatialData:
        # subsetting the data for debug/development purposes
        min_coordinate = [12790, 12194]
        max_coordinate = [15100, 15221]
        merged = merged.query.bounding_box(
            min_coordinate=min_coordinate,
            max_coordinate=max_coordinate,
            axes=["y", "x"],
            target_coordinate_system="aligned",
        )
        if False:
            from napari_spatialdata import Interactive

            Interactive(merged)
            os._exit(0)
        return merged

    @staticmethod
    def _napari_visualization(merged: SpatialData, dataset: ImageTilesDataset):
        ##
        # visualization (only if using the small dataset, otherwise the code will be too slow)
        tiles = []
        for i, tile in enumerate(tqdm(dataset)):
            if i > 100:
                break
            tiles.append(tile)

        tiles_sdata = SpatialData(
            images={
                "full": merged.images["morphology_mip"]
                # "full": merged.images["CytAssist_FFPE_Human_Breast_Cancer_full_image"]
            }
            | {f"tile_{j}": tile for j, (tile, _) in enumerate(tiles)}
        )

        Interactive([tiles_sdata, merged])

    def __len__(self):
        return len(self.dataset0)

    def __getitem__(self, idx):
        tile0, row = self.dataset0[idx]
        tile1, _ = self.dataset1[idx]
        expected_category = row.obs[self.CATEGORICAL_COLUMN].values[0]
        expected_category = self.categories.index(expected_category)
        clonal_one_hot = F.one_hot(
            torch.tensor(expected_category), num_classes=len(self.categories)
        )
        tile0_numpy = tile0.data.compute()
        tile0_numpy = (tile0_numpy - self.visium_mean[:, None, None]) / self.visium_std[:, None, None]
        tile0_numpy = tile0_numpy.astype(np.float32)
        tile1_numpy = tile1.data.compute()
        tile1_numpy = (tile1_numpy - self.xenium_mean) / self.xenium_std
        tile1_numpy = tile1_numpy.astype(np.float32)
        stacked_tile = np.vstack([tile0_numpy, tile1_numpy])
        if self.transform is not None:
            stacked_tile = self.transform(stacked_tile)
        return stacked_tile, clonal_one_hot.float()


class TilesDataModule(LightningDataModule):
    def __init__(
        self,
        batch_size: int,
        num_workers: int,
        train_transform: torchvision.transforms.Compose,
    ):
        super().__init__()
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.train_transform = train_transform
        self.dataset = None
        self.train = None
        self.val = None
        self.test = None
        self.num_classes = None
        self.in_channels = None
        self.dataloader_kwargs = {
            "pin_memory": True if torch.cuda.is_available() else False,
            "persistent_workers": True if torch.cuda.is_available() else False,
        }

    def setup(self, stage=None):
        self.dataset = TileDataset()
        n_train = int(len(self.dataset) * 0.7)
        n_val = int(len(self.dataset) * 0.2)
        n_test = len(self.dataset) - n_train - n_val
        self.train, self.val, self.test = torch.utils.data.random_split(
            self.dataset,
            [n_train, n_val, n_test],
            generator=torch.Generator().manual_seed(42),
        )
        self.train.dataset.transform = self.train_transform
        self.train.dataset.transform = self.train_transform
        self.num_classes = len(self.dataset.categories)
        self.in_channels = self.dataset[0][0].shape[0]

    def train_dataloader(self):
        return DataLoader(
            self.train,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            shuffle=True,
            **self.dataloader_kwargs,
        )

    def val_dataloader(self):
        return DataLoader(
            self.val,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            shuffle=False,
            **self.dataloader_kwargs,
        )

    def test_dataloader(self):
        return DataLoader(
            self.test,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            shuffle=False,
            **self.dataloader_kwargs,
        )

    def predict_dataloader(self):
        return DataLoader(
            self.dataset,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            shuffle=False,
            **self.dataloader_kwargs,
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
