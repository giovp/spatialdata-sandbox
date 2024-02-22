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
            "xenium_1.0.2_io/to_zarr.py",
            "xenium_1.3.0_io/to_zarr.py",
            "xenium_1.4.0_io/to_zarr.py",
            "xenium_1.5.0_io/to_zarr.py",
            "xenium_1.6.0_io/to_zarr.py",
            "xenium_1.7.0_io/to_zarr.py",
            "xenium_2.0.0_io/to_zarr.py",
        ],
    )
    def test_xenium_to_zarr_scripts(self, script):
        self.execute_script(path / script)
