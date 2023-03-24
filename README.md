# Foldcomp

<p align="center">
<img src="https://raw.githubusercontent.com/steineggerlab/foldcomp/master/.github/img/foldcomp_strong_marv.png" max-height="300px" height="300" display="block" margin-left="auto" margin-right="auto" display="block"/>
</p>
Foldcomp compresses protein structures with torsion angles effectively. It compresses the backbone atoms to 8 bytes and the side chain to additionally 4-5 byes per residue, thus an averaged-sized protein of 350 residues requires ~6kb.

Foldcomp efficient compressed format stores protein structures requiring only 13 bytes per residue, which reduces the required storage space by an order of magnitude compared to saving 3D coordinates directly. We achieve this reduction by encoding the torsion angles of the backbone as well as the side-chain angles in a compact binary file format (FCZ).

> Foldcomp currently only supports compression of single chain PDB files
<br clear="right"/>

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/steineggerlab/foldcomp/master/.github/img/format_benchmark_dark.png">
  <img src="https://raw.githubusercontent.com/steineggerlab/foldcomp/master/.github/img/format_benchmark_light.png" alt="Left panel: Foldcomp data format, saving amino acid residue in 13 byte. Top right panel:  Foldcomp decompression is as fast as gzip. Bottom right panel: Foldcomp compression ratio is higher than pulchra and gzip." max-width="720px" max-height="400px" width="auto" height="auto">
</picture>
</p>

## Publications

[Hyunbin Kim, Milot Mirdita, and Martin Steinegger. Foldcomp: a library and format for compressing and indexing large protein structure sets. bioRxiv, doi:10.1101/2022.12.09.519715  (2022)](https://www.biorxiv.org/content/10.1101/2022.12.09.519715v1)

## Usage

### Installing Foldcomp

```
# Install Foldcomp Python package
pip install foldcomp

# Download static binaries for Linux
wget https://mmseqs.com/foldcomp/foldcomp-linux-x86_64.tar.gz

# Download static binaries for Linux (ARM64)
wget https://mmseqs.com/foldcomp/foldcomp-linux-arm64.tar.gz

# Download binary for macOS
wget https://mmseqs.com/foldcomp/foldcomp-macos-universal.tar.gz

# Download binary for Windows (x64)
wget https://mmseqs.com/foldcomp/foldcomp-windows-x64.zip
```

### Executable
```
# Compression
foldcomp compress <pdb_file|cif_file> [<fcz_file>]
foldcomp compress [-t number] <pdb_dir|cif_dir> [<fcz_dir>]

# Decompression
foldcomp decompress <fcz_file> [<pdb_file>]
foldcomp decompress [-t number] <fcz_dir> [<pdb_dir>]

# Extraction of sequence or pLDDT
foldcomp extract [--plddt|--fasta] <fcz_file> [<txt_file|fasta_file>]
foldcomp extract [--plddt|--fasta] [-t number] <fcz_dir|tar> [<output_dir>]

# Check
foldcomp check <fcz_file>
foldcomp check [-t number] <fcz_dir|tar>

# RMSD
foldcomp rmsd <pdb1|cif1> <pdb2|cif2>

# Options
 -h, --help           print this help message
 -t, --threads        threads for (de)compression of folders/tar files [default=1]
 -r, --recursive      recursively look for files in directory [default=0]
 -f, --file           input is a list of files [default=0]
 -a, --alt            use alternative atom order [default=false]
 -b, --break          interval size to save absolute atom coordinates [default=25]
 -z, --tar            save as tar file [default=false]
 -d, --db             save as database [default=false]
 -y, --overwrite          overwrite existing files [default=false]
 --skip-discontinuous skip PDB with with discontinuous residues (only batch compression)
 --plddt              extract pLDDT score (only for extraction mode)
 --fasta              extract amino acid sequence (only for extraction mode)
 --no-merge           do not merge output files (only for extraction mode)
 --time               measure time for compression/decompression
```

### Downloading Databases
We offer prebuilt databases for multiple large sets of predicted protein structures and a Python helper to download the database files.

You can download the AlphaFoldDB Swiss-Prot with the following command:
```
python -c "import foldcomp; foldcomp.setup('afdb_swissprot_v4');
```

Currently we offer the following databases:
* [ESMAtlas v2023_02](https://esmatlas.com/): `foldcomp.setup('esmatlas_v2023_02')`
* ESMAtlas high-quality: `foldcomp.setup('highquality_clust30')`

  **Note:** We skipped all structures with discontinous residues or other issues.
   Here is a list with the affected predictions;
   [high-quality](https://foldcomp.steineggerlab.workers.dev/highquality_clust30_issues.txt) (~100k),
   [v2023_02](https://foldcomp.steineggerlab.workers.dev/esmatlas_v2023_02.err.log) (~10k)

* [AlphaFoldDB Uniprot](https://alphafold.ebi.ac.uk/): `foldcomp.setup('afdb_uniprot_v4')`
* AlphaFoldDB Swiss-Prot: `foldcomp.setup('afdb_swissprot_v4')`
* [AlphaFoldDB Cluster Representatives](https://afdb-cluster.steineggerlab.workers.dev/): `foldcomp.setup('afdb_rep_v4')`
* AlphaFoldDB Cluster Representatives (Dark Clusters): `foldcomp.setup('afdb_rep_dark_v4')`

If you want other prebuilt datasets, please get in touch with us through our [GitHub issues](https://github.com/steineggerlab/foldcomp/issues).

If you have issues downloading the databases you can navigate directly to our [download server](https://foldcomp.steineggerlab.workers.dev/) and download the required files. E.g. `afdb_uniprot_v4`, `afdb_uniprot_v4.index`, `afdb_uniprot_v4.dbtype`, `afdb_uniprot_v4.lookup`, and optionally `afdb_uniprot_v4.source`.

### Python API

You can find more in-depth examples of using Foldcomp's Python interface in the example notebook:
<a href="https://colab.research.google.com/github/steineggerlab/foldcomp/blob/master/foldcomp-py-examples.ipynb" target="_blank" rel="noopener"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

```py
import foldcomp
# 01. Handling a FCZ file
# Open a fcz file
with open("test/compressed.fcz", "rb") as fcz:
  fcz_binary = fcz.read()

  # Decompress
  (name, pdb) = foldcomp.decompress(fcz_binary) # pdb_out[0]: file name, pdb_out[1]: pdb binary string

  # Save to a pdb file
  with open(name, "w") as pdb_file:
    pdb_file.write(pdb)

  # Get data as dictionary (v0.0.3)
  data_dict = foldcomp.get_data(fcz_binary) # foldcomp.get_data(pdb) also works
  # Keys: phi, psi, omega, torsion_angles, residues, bond_angles, coordinates
  data_dict["phi"] # phi angles (C-N-CA-C)
  data_dict["psi"] # psi angles (N-CA-C-N)
  data_dict["omega"] # omega angles (CA-C-N-CA)
  data_dict["torsion_angles"] # torsion angles of the backbone as list (phi + psi + omega)
  data_dict["bond_angles"] # bond angles of the backbone as list
  data_dict["residues"] # amino acid residues as string
  data_dict["coordinates"] # coordinates of the backbone as list

# 02. Iterate over a database of FCZ files
# Open a foldcomp database
ids = ["d1asha_", "d1it2a_"]
with foldcomp.open("test/example_db", ids=ids) as db:
  # Iterate through database
  for (name, pdb) in db:
      # save entries as seperate pdb files
      with open(name + ".pdb", "w") as pdb_file:
        pdb_file.write(pdb)
```

## Subsetting Databases
If you are dealing with millions of entries, we recommend using `createsubdb` command
of [mmseqs2](https://mmseqs.com) to subset databases.
The following commands can be used to subset the AlphaFold Uniprot DB with given IDs.
```sh
# mmseqs createsubdb --subdb-mode 0 --id-mode 1 id_list.txt input_foldcomp_db output_foldcomp_db
mmseqs createsubdb --subdb-mode 0 --id-mode 1 id_list.txt afdb_uniprot_v4 afdb_subset
```
Please note that the IDs in afdb_uniprot_v4 are in the format `AF-A0A5S3Y9Q7-F1-model_v4` .

## Community Contributions
* [PyMOL Plugin for reading Foldcomp files](https://github.com/yakomaxa/load_fcz_PyMOL) by @yakomaxa

## Contributor
<a href="https://github.com/steineggerlab/foldcomp/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/foldcomp" />
</a>

