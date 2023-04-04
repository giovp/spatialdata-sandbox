import os

import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor
from pytorch_lightning.callbacks.progress import TQDMProgressBar
from pytorch_lightning.loggers import CSVLogger

import torchvision
import torch

from examples.xenium_and_visium_data import TilesDataModule, DenseNetModel
import numpy as np


def random_horizontal_flip(image: np.ndarray, p: float = 0.5) -> np.ndarray:
    """
    Randomly apply horizontal flip to a cyx numpy array image.

    :param image: cyx numpy array image (c = channels, y = height, x = width)
    :param p: probability of applying the horizontal flip (default: 0.5)
    :return: numpy array of the horizontally flipped image
    """
    if np.random.rand() < p:
        flipped_image = np.flip(image, axis=-1).copy()  # Flip along the x-axis
        return flipped_image
    else:
        return image


if __name__ == "__main__":
    pl.seed_everything(7)

    PATH_DATASETS = os.environ.get("PATH_DATASETS", ".")
    BATCH_SIZE = 4096 if torch.cuda.is_available() else 2
    NUM_WORKERS = 10 if torch.cuda.is_available() else 6
    print(f"Using {BATCH_SIZE} batch size.")
    print(f"Using {NUM_WORKERS} workers.")

    train_transform = torchvision.transforms.Compose([random_horizontal_flip])

    tiles_data_module = TilesDataModule(
        batch_size=BATCH_SIZE,
        num_workers=NUM_WORKERS,
        train_transform=train_transform,
    )

    tiles_data_module.setup()
    train_dl = tiles_data_module.train_dataloader()
    val_dl = tiles_data_module.val_dataloader()
    test_dl = tiles_data_module.test_dataloader()

    model = DenseNetModel(
        learning_rate=1e-5,
        in_channels=tiles_data_module.in_channels,
        num_classes=tiles_data_module.num_classes,
    )
    import logging

    logging.basicConfig(level=logging.INFO)

    trainer = pl.Trainer(
        max_epochs=5,
        accelerator="auto",
        devices=1 if torch.cuda.is_available() else None,  # limiting got iPython runs
        logger=CSVLogger(save_dir="logs/"),
        callbacks=[
            LearningRateMonitor(logging_interval="step"),
            TQDMProgressBar(refresh_rate=5),
        ],
        log_every_n_steps=50,
    )

    trainer.fit(model, datamodule=tiles_data_module)
    trainer.test(model, datamodule=tiles_data_module)
