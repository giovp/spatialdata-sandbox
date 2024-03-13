import pathlib


path = pathlib.Path(__file__).parent.parent.resolve()

import subprocess
import pytest


class TestToZarrScripts:
    def execute_script(self, script):
        result = subprocess.run(
            ["python", script],
            capture_output=True,
            cwd=script.parent,
            text=True,
        )

        # Check that script finished successfully (returncode 0)
        if result.returncode != 0:
            print(
                f"Script {script} failed with output:\n{result.stdout}\n{result.stderr}"
            )
        assert result.returncode == 0

    @pytest.mark.parametrize(
        "script",
        [
            "visium_1.0.0_io/to_zarr.py",
            "visium_1.1.0_io/to_zarr.py",
            "visium_1.2.0_io/to_zarr.py",
            "visium_1.3.0_io/to_zarr.py",
            # "visium_1.3.1_io/to_zarr.py",  # not supported due to the format appearing very different on disk
            "visium_2.0.0_1_io/to_zarr.py",
            "visium_2.0.0_2_io/to_zarr.py",
            "visium_2.0.1_io/to_zarr.py",
            "visium_2.1.0_1_io/to_zarr.py",
            "visium_2.1.0_2_io/to_zarr.py",
            "visium_2.1.0_3_io/to_zarr.py",
        ],
    )
    def test_xenium_to_zarr_scripts(self, script):
        self.execute_script(path / script)
