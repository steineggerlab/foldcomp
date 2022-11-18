from skbuild import setup

setup(
    name="foldcomp",
    version="0.0.3",
    description="Foldcomp compresses protein structures with torsion angles effectively. It compresses the backbone atoms to 8 bytes and the side chain to additionally 4-5 byes per residue, an averaged-sized protein of 350 residues requires ~4.2kb. Foldcomp is a C++ library with Python bindings.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Milot Mirdita <milot@mirdita.de>, Hyunbin Kim <khb7840@gmail.com>, Martin Steinegger <themartinsteinegger@gmail.com>",
    license="GPLv3",
    cmake_args=["-DBUILD_PYTHON:BOOL=ON"],
    python_requires=">=3.7",
    packages=["foldcomp"],
    include_package_data=False,
    install_requires=[
        "httpx >= 0.23.0",
    ],
    extras_require={"test": ["pytest"]},
)
