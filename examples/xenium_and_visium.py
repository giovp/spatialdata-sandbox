##
import os
from anndata import AnnData
import torch.nn.functional as F
import numpy as np
from spatial_image import SpatialImage
import spatialdata as sd
from napari_spatialdata import Interactive
from spatialdata import SpatialData
from spatialdata.transformations import get_transformation, set_transformation
from spatialdata.dataloader.datasets import ImageTilesDataset
from tqdm import tqdm
from spatialdata import transform
from typing import Dict, Callable
import os
import scanpy as sc

from pl_bolts.datamodules import CIFAR10DataModule
from pl_bolts.transforms.dataset_normalizations import cifar10_normalization
from monai.networks.nets import DenseNet121
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor
from pytorch_lightning.callbacks.progress import TQDMProgressBar
from pytorch_lightning.loggers import CSVLogger
import torch
from torch.nn import CrossEntropyLoss
from torch.optim import Adam
import torchvision
import torch
import torchvision.transforms as transforms
from torchvision.datasets import ImageFolder
from torch.utils.data import DataLoader
from pytorch_lightning import LightningDataModule

## TODO:
# replace clonal information with one-hot encoding
# setup tensorboard to train model
# train model on small data
# deploy on GPU machine
# train model
# hyperparamers optimization?
# plot result with napari

##
SMALL_DATASET = False

# loading the Xenium and Visium SpatialData objects
print("current working directory:", os.getcwd())
SPATIALDATA_SANDBOX_PATH = "."
GENERATED_DATA_PATH = os.path.join(
    SPATIALDATA_SANDBOX_PATH, "generated_data/xenium_visium_integration"
)
assert os.path.isdir(
    GENERATED_DATA_PATH
), f"{GENERATED_DATA_PATH} not found, please use symlinks to make it available"

XENIUM_SDATA_PATH = os.path.join(SPATIALDATA_SANDBOX_PATH, "xenium_rep1_io/data.zarr")
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
# 1. ensuring the data is aligned (i.e. that the coordinate system "aligned" is present both for the Xenium and Visium
#    data
# 2. removing unnecessary coordinate systems
for sdata in [xenium_sdata, visium_sdata]:
    aligned = False
    all_elements = list(sdata._gen_elements())
    for element_type, name, element in all_elements:
        all_transformations = get_transformation(element, get_all=True)
        if "aligned" in all_transformations:
            set_transformation(
                element, {"aligned": all_transformations["aligned"]}, set_all=True
            )
            aligned = True
        else:
            print(
                "element",
                element_type,
                name,
                'is not needed because it has no transformation to the "aligned" coordinate system, removing it',
            )
            del sdata.__getattribute__(element_type)[name]
    if not aligned:
        raise ValueError(
            "No transformation to the 'aligned' coordinate system was found in any element of the spatial data. "
            "Please run the appropriate notebook to aligned the data before running the current code."
        )

##
# create a dataset with tiles, one for each Visium circle, but around the Xenium image
merged = sd.SpatialData(
    images={
        "xenium": xenium_sdata.images["morphology_mip"],
        "visium": visium_sdata.images["CytAssist_FFPE_Human_Breast_Cancer_full_image"],
    },
    shapes={
        "CytAssist_FFPE_Human_Breast_Cancer": visium_sdata.shapes[
            "CytAssist_FFPE_Human_Breast_Cancer"
        ]
    },
    table=visium_sdata.table,
)

s = adata.obs["clone"]
s.index = sdata.table.obs.index
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
        tile, region_name, region_index = self.tile_dataset[idx]
        clonal_information = (
            self.tile_dataset.sdata.table[region_index].obs["clone"].values[0]
        )
        clonal_information = self.categories.index(clonal_information)
        clonal_one_hot = F.one_hot(
            torch.tensor(clonal_information), num_classes=len(self.categories)
        )
        return tile, clonal_one_hot.float()


clonal_tile_dataset = ClonalTileDataset(dataset)

##
# visualization (only if using the small dataset, otherwise the code will be too slow)
if SMALL_DATASET:
    tiles = []
    for i, tile in enumerate(tqdm(dataset)):
        tiles.append(tile)

    tiles_sdata = SpatialData(
        images={
            "full": visium_sdata.images["CytAssist_FFPE_Human_Breast_Cancer_full_image"]
        }
        | {
            f"tile_{region_name}_{region_index}": tile
            for (tile, region_name, region_index) in tiles
        }
    )

    Interactive([tiles_sdata, merged])
##

# set up a classifier on the tiles
pl.seed_everything(7)

PATH_DATASETS = os.environ.get("PATH_DATASETS", ".")
BATCH_SIZE = 256 if torch.cuda.is_available() else 64
# NUM_WORKERS = int(os.cpu_count() / 2)
NUM_WORKERS = 1


##
def spatial_image_to_tensor(image: SpatialImage):
    return torch.from_numpy(image.data.compute().astype(np.float16))


train_transform = torchvision.transforms.Compose(
    [
        # torchvision.transforms.RandomCrop(32, padding=4),
        spatial_image_to_tensor,
        torchvision.transforms.RandomHorizontalFlip(),
    ]
)

test_transform = torchvision.transforms.Compose(
    [
        spatial_image_to_tensor,
    ]
)


class TilesDataModule(LightningDataModule):
    def __init__(
        self,
        batch_size: int,
        num_workers: int,
        train_transform: torchvision.transforms.Compose,
        test_transform: torchvision.transforms.Compose,
        dataset: torch.utils.data.Dataset,
    ):
        super().__init__()
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.train_transform = train_transform
        self.test_transform = test_transform
        self.dataset = dataset
        self.train = None
        self.val = None
        self.test = None

    def setup(self, stage=None):
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


tiles_data_module = TilesDataModule(
    batch_size=BATCH_SIZE,
    num_workers=NUM_WORKERS,
    train_transform=train_transform,
    test_transform=test_transform,
    dataset=clonal_tile_dataset,
)

tiles_data_module.setup()
train_dl = tiles_data_module.train_dataloader()
val_dl = tiles_data_module.val_dataloader()
test_dl = tiles_data_module.test_dataloader()


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
        return self._compute_loss_from_batch(batch=batch, batch_idx=batch_idx)

    def test_step(self, batch, batch_idx):
        imgs, labels = batch
        preds = self.model(imgs).argmax(dim=-1)
        acc = (labels == preds).float().mean()
        # By default logs it per epoch (weighted average over batches), and returns it afterwards
        self.log("test_acc", acc)

    def configure_optimizers(self) -> Adam:
        return Adam(self.model.parameters(), lr=self.hparams.learning_rate)


if __name__ == "__main__":
    t = train_dl.__iter__().__next__()
    print(t)
    model = DenseNetModel(
        learning_rate=0.05, in_channels=1, num_classes=len(_get_clonal_classes(adata))
    )

    trainer = pl.Trainer(
        max_epochs=30,
        accelerator="auto",
        devices=1 if torch.cuda.is_available() else None,  # limiting got iPython runs
        logger=CSVLogger(save_dir="logs/"),
        callbacks=[
            LearningRateMonitor(logging_interval="step"),
            TQDMProgressBar(refresh_rate=10),
        ],
    )

    trainer.fit(model, datamodule=tiles_data_module)
    trainer.test(model, datamodule=tiles_data_module)
