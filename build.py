from skbuild import setup, constants
from pathlib import Path
import shutil
import sys

# set debug mode for scikit-build
# sys.argv[2:2] = ["--build-type", "Debug"]

setup(cmake_args=["-DBUILD_PYTHON:BOOL=ON"], script_args=["build"])
src_dir = Path(constants.CMAKE_INSTALL_DIR()) / "foldcomp"
dest_dir = Path("foldcomp")


def remove_files(target_dir: Path, pattern: str) -> None:
    """Delete files matched with a glob pattern in a directory tree."""
    for path in target_dir.glob(pattern):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def copy_files(src_dir: Path, dest_dir: Path, pattern: str) -> None:
    """Copy files matched with a glob pattern in a directory tree to another."""
    for src in src_dir.glob(pattern):
        dest = dest_dir / src.relative_to(src_dir)
        if src.is_dir():
            # NOTE: inefficient if subdirectories also match to the pattern.
            copy_files(src, dest, "*")
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)


# Delete C-extensions copied in previous runs, just in case.
remove_files(dest_dir, "**/*.pyd")
remove_files(dest_dir, "**/*.so")

# Copy built C-extensions back to the project.
copy_files(src_dir, dest_dir, "**/*.pyd")
copy_files(src_dir, dest_dir, "**/*.so")
