import os
import subprocess
import urllib.request
import zipfile

from tqdm import tqdm


class TqdmDownload(tqdm):
    def __init__(self, *args, **kwargs):
        kwargs = dict(kwargs)
        kwargs.update({"unit": "B", "unit_scale": True, "unit_divisor": 1024})
        super().__init__(*args, **kwargs)

    def update_to(self, nblocks=1, blocksize=1, total=-1):
        self.total = total
        self.update(nblocks * blocksize - self.n)


def is_aria2c_installed():
    rc = subprocess.call(["which", "aria2c"])
    if rc == 0:
        return True
    else:
        return False


def download(url, outfile, desc):
    # aria2c (maybe) is faster than urllib
    if not is_aria2c_installed():
        with TqdmDownload(desc="downloading " + desc) as t:
            urllib.request.urlretrieve(url, outfile, t.update_to)
    else:
        subprocess.check_call(
            [
                "aria2c",
                "-x",
                "4",
                "--dir",
                os.path.dirname(outfile),
                "-o",
                os.path.basename(outfile),
                url,
            ]
        )


def unzip(file, outdir=None, files=None, rm=True):
    if outdir is None:
        outdir = os.getcwd()
    else:
        os.makedirs(outdir, exist_ok=True)
    zfile = zipfile.ZipFile(file)
    if files is not None:
        for f in files:
            zfile.extract(f, outdir)
    else:
        zfile.extractall(outdir)
    zfile.close()
    if rm:
        os.unlink(file)
