# Run from the repo root:
#   asv run --config xenium_2.0.0_io/asv.conf.json HEAD^!
#
# IMPORTANT: use HEAD^! (not HEAD). HEAD^! means "this commit only".
# `asv run HEAD` iterates ALL commits reachable from HEAD (~188 commits),
# running hundreds of benchmarks instead of the intended 2.
#
# To also benchmark the locally installed spatialdata (e.g. a dev branch),
# run a second pass using the existing environment (no commit range allowed):
#   asv run --config xenium_2.0.0_io/asv.conf.json --environment existing:python3.13
#
# Show results in the terminal (all environments side by side):
#   asv show --config xenium_2.0.0_io/asv.conf.json HEAD
#
# Build and serve an interactive HTML report:
#   asv publish --config xenium_2.0.0_io/asv.conf.json
#   asv preview --config xenium_2.0.0_io/asv.conf.json
#   # then open http://127.0.0.1:8080 in your browser

import tempfile
from pathlib import Path

import spatialdata_io

# Data is at xenium_2.0.0_io/data/, this file is at benchmarks/
DATA_PATH = Path(__file__).parent.parent / "xenium_2.0.0_io" / "data"


class XeniumWriteSuite:
    """
    Benchmark sdata.write() for the xenium 2.0.0 dataset.

    Compares spatialdata v0.7.2 vs v0.7.3a0 (configured via asv.conf.json matrix).
    Only the write step is timed; xenium parsing runs in setup() and is excluded.
    """

    # Write is slow: run once per benchmark session
    number = 1
    repeat = 1
    warmup_time = 0
    timeout = 600

    def setup(self):
        if not DATA_PATH.exists():
            raise FileNotFoundError(
                f"Xenium data not found at {DATA_PATH}. "
                "Run download.py first."
            )
        self.sdata = spatialdata_io.xenium(
            cells_boundaries = False,
            nucleus_boundaries = False,
            cells_as_circles = False,
            cells_labels = False,
            nucleus_labels = False,
            transcripts = False,
            morphology_mip = False,
            morphology_focus = True,
            aligned_images = False,
            cells_table = False,
            gex_only = False,
            path=str(DATA_PATH),
        )
        self._tmpdir = tempfile.TemporaryDirectory(prefix="asv_xenium_write_")

    def teardown(self):
        self._tmpdir.cleanup()

    def time_write(self):
        self.sdata.write(Path(self._tmpdir.name) / "data.zarr")

    def peakmem_write(self):
        self.sdata.write(Path(self._tmpdir.name) / "data.zarr")
