import pytorch_lightning as pl
from pytorch_lightning.callbacks.progress import TQDMProgressBar
import torchvision
from examples.xenium_and_visium0 import random_horizontal_flip
import torch
import os
import re
import numpy as np
import pandas as pd

from examples.xenium_and_visium_data import (
    TilesDataModule,
    DenseNetModel,
)


def find_latest_checkpoint(logs_dir: str) -> tuple[str, str]:
    assert os.path.isdir(logs_dir), logs_dir
    highest_version = -1
    highest_version_ckpt = None

    # Iterate through directories in the lightning_logs folder
    for entry in os.scandir(logs_dir):
        if entry.is_dir() and entry.name.startswith("version_"):
            # Extract the version number using a regular expression
            match = re.match(r"version_(\d+)", entry.name)
            if match:
                version_number = int(match.group(1))

                # Update the highest version number and corresponding checkpoint path
                if version_number > highest_version:
                    highest_version = version_number
                    checkpoints_folder = os.path.join(entry.path, "checkpoints")
                    valid = True
                    if os.path.exists(checkpoints_folder):
                        ckpt_files = [
                            f
                            for f in os.listdir(checkpoints_folder)
                            if f.endswith(".ckpt")
                        ]
                        if len(ckpt_files) == 1:
                            highest_version_ckpt = os.path.join(
                                entry.path, "checkpoints", ckpt_files[0]
                            )
                        else:
                            valid = False
                    else:
                        valid = False
                    if not valid:
                        highest_version_ckpt = None

    if highest_version == -1 or highest_version_ckpt is None:
        raise ValueError(
            "No valid version directory found in the specified logs directory."
        )

    assert os.path.isfile(highest_version_ckpt), highest_version_ckpt
    return highest_version_ckpt


if __name__ == "__main__":
    print(os.getcwd())
    latest_checkpoint = find_latest_checkpoint(logs_dir="logs/lightning_logs/")
    # latest_checkpoint = find_latest_checkpoint(logs_dir="/home/l989o/pycharm_project_143/logs/lightning_logs/")
    print("loading checkpoint:", latest_checkpoint)

    model = DenseNetModel.load_from_checkpoint(latest_checkpoint)

    # disable randomness, dropout, etc...
    model.eval()

    BATCH_SIZE = 4096 if torch.cuda.is_available() else 64
    NUM_WORKERS = 10 if torch.cuda.is_available() else 6
    #
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

    trainer = pl.Trainer(
        accelerator="auto",
        devices=1 if torch.cuda.is_available() else None,  # limiting got iPython runs
        callbacks=[
            TQDMProgressBar(refresh_rate=10),
        ],
    )

    predictions = trainer.predict(datamodule=tiles_data_module, model=model)
    predictions = torch.cat(predictions, dim=0)

    print(np.unique(predictions.detach().cpu().numpy(), return_counts=True))

    ##
    p = predictions.detach().cpu().numpy()
    categories = tiles_data_module.dataset.categories
    predicted_celltype_major = []
    for i in p:
        predicted_celltype_major.append(categories[i])
    s = pd.Series(predicted_celltype_major)
    categorical = pd.Categorical(s, categories=categories)

    sdata = tiles_data_module.dataset.dataset0.sdata
    categorical.index = sdata.table.obs.index
    sdata.table.obs["predicted_celltype_major"] = categorical
    ##

    from napari_spatialdata import Interactive

    Interactive(sdata)
