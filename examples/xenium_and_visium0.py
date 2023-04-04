import os

import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor
from pytorch_lightning.callbacks.progress import TQDMProgressBar
from pytorch_lightning.loggers import CSVLogger

import torchvision
import torch

from examples.xenium_and_visium_data import TilesDataModule, DenseNetModel, _get_clonal_classes

if __name__ == "__main__":
    pl.seed_everything(7)

    PATH_DATASETS = os.environ.get("PATH_DATASETS", ".")
    BATCH_SIZE = 512 if torch.cuda.is_available() else 64
    NUM_WORKERS = 64

    train_transform = torchvision.transforms.Compose(
        [
            torchvision.transforms.RandomHorizontalFlip(),
        ]
    )

    test_transform = torchvision.transforms.Compose(
        [
        ]
    )
    tiles_data_module = TilesDataModule(
        batch_size=BATCH_SIZE,
        num_workers=NUM_WORKERS,
        train_transform=train_transform,
        test_transform=test_transform,
    )

    tiles_data_module.setup()
    train_dl = tiles_data_module.train_dataloader()
    val_dl = tiles_data_module.val_dataloader()
    test_dl = tiles_data_module.test_dataloader()

    adata = tiles_data_module.dataset.tile_dataset.sdata.table
    model = DenseNetModel(
        learning_rate=0.05, in_channels=1, num_classes=len(_get_clonal_classes(adata))
    )

    trainer = pl.Trainer(
        max_epochs=5,
        accelerator="auto",
        devices=1 if torch.cuda.is_available() else None,  # limiting got iPython runs
        logger=CSVLogger(save_dir="logs/"),
        callbacks=[
            LearningRateMonitor(logging_interval="step"),
            TQDMProgressBar(refresh_rate=10),
        ],
        log_every_n_steps=7,
    )

    trainer.fit(model, datamodule=tiles_data_module)
    trainer.test(model, datamodule=tiles_data_module)
