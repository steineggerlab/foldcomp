from skbuild import setup

setup(
    name="foldcomp",
    version="0.0.1",
    description="Foldcomp compresses protein structures with torsion angles effectively. It compresses the backbone atoms to 8 bytes and the side chain to additionally 4-5 byes per residue, an averaged-sized protein of 350 residues requires ~4.2kb. Foldcomp is a C++ library with Python bindings.",
    author="Milot Mirdita <milot@mirdita.de>",
    license="MIT",
    cmake_args=["-DBUILD_PYTHON:BOOL=ON"],
    python_requires=">=3.7",
    packages=["foldcomp"],
    include_package_data=False,
    install_requires=[
        "aiohttp >= 3.8.3",
    ],
    extras_require={"test": ["pytest"]},
)
