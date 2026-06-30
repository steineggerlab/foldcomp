from pathlib import Path
import shutil

from skbuild import setup


def _remove_stale_skbuild_cache() -> None:
    """Remove cached scikit-build build trees before configuring CMake.

    This repo is often built from a workspace that already contains generated
    `_skbuild` artifacts. If the cached CMake build tree was created with a
    different generator, CMake refuses to reconfigure it. Clearing the cached
    build tree keeps `pip`/`uv` installs reproducible from a dirty checkout.
    """

    project_root = Path(__file__).resolve().parent
    skbuild_root = project_root / "_skbuild"
    if not skbuild_root.exists():
        return

    for build_dir in skbuild_root.glob("*/cmake-build"):
        shutil.rmtree(build_dir, ignore_errors=True)


_remove_stale_skbuild_cache()

setup(
    name="foldcomp",
    version="1.0.0",
    description="Foldcomp compresses protein structures with torsion angles effectively. It compresses the backbone atoms to 8 bytes and the side chain to additionally 4-5 byes per residue, an averaged-sized protein of 350 residues requires ~4.2kb. Foldcomp is a C++ library with Python bindings.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Milot Mirdita <milot@mirdita.de>, Hyunbin Kim <khb7840@gmail.com>, Martin Steinegger <themartinsteinegger@gmail.com>",
    license="MIT",
    cmake_args=["-DBUILD_PYTHON:BOOL=ON"],
    python_requires=">=3.7",
    packages=["foldcomp"],
    include_package_data=False,
    install_requires=[
        "httpx >= 0.23.0",
    ],
    extras_require={"test": ["pytest"]},
)
