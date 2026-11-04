import logging
import pathlib
import subprocess

logger: logging.Logger = logging.getLogger(__name__)


class GzipManager:
    """Context manager for creating an unzipped file if needed."""

    def __init__(self, pth: str):
        self.pth = pathlib.Path(pth)
        # Determine if this is gzipped, based on suffix
        self.gzipped = self.pth.suffix == ".gz"
        self.pth_unzipped = self.pth.with_suffix("")

    def __enter__(self):
        if self.gzipped:
            logger.debug(f"Unzipping file {self.pth}")
            cmd = f"gunzip -c {self.pth} > {self.pth_unzipped}"
            res = subprocess.run(cmd, shell=True, check=True)
            if res.returncode != 0:
                raise RuntimeError(
                    f"gunzip failed with return code {res.returncode}"
                )
            return pathlib.Path(f"{self.pth_unzipped}")
        else:
            logger.debug(f"File {self.pth} is already unzipped")
            return self.pth

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.gzipped:
            logger.debug(f"Cleaning up decompressed file {self.pth_unzipped}")
            self.pth_unzipped.unlink()
