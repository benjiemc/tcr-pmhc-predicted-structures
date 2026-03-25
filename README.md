# TCR:pMHC Predicted Structures

This project creates synthetic TCR:pMHC structures from sequence data using state of the art deep learning models.

## Downloading required data

This project requires several datasets before it can be run. These include:

* The AlphaFold 3 weights: these can be downloaded following the instructions [here](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#obtaining-model-parameters).
* The AlphaFold 3 genetic databases: these can be downloaded following the instructions [here](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#obtaining-genetic-databases).
* The ColabFold databases: these can be downloaded following the instructions [here](https://github.com/sokrypton/colabfold?tab=readme-ov-file#generating-msas-for-large-scale-structurecomplex-predictions).
* The TCRmodel2 databases: these should be downloaded following the instructions [here](https://github.com/piercelab/tcrmodel2?tab=readme-ov-file#download-database) and placed in a new directory.

Once complete, the environment variables `ALPHAFOLD3_WEIGHTS_DIR`, `ALPHAFOLD_DATA_DIR`, `COLABFOLD_DATA_DIR`, and `TCRMODEL2_ALPHAFOLD_DATA_DIR` should be set to the corresponding locations.

## Installation and Setup

To install this project, first clone the repo:

```
git clone git@github.com:benjiemc/tcr-pmhc-predicted-structures.git
cd tcr-pmhc-predicted-structures/
```

The following envrionment variables need to be set during the initialisation stage:

- `ALPHAFOLD3_WEIGHTS_DIR`: Path to the AF3 weights directory
- `ALPHAFOLD_DATA_DIR`: Path to the AlphaFold3 data directory
- `COLABFOLD_DATA_DIR`: Path to the ColabFold data directory
- `TCRMODEL2_ALPHAFOLD_DATA_DIR`: Path to the TCRModel2 data directory

> ### Using MODELLER to fix missing residues and atoms
>
> Parts of the data pipeline for this project use [MODELLER](https://salilab.org/modeller/) to fix structures missing atoms and residues. A license is required to use MODELLER. Please [register for a license](https://salilab.org/modeller/registration.html) and then set the license key as an evironment variable using the following command: `export KEY_MODELLER='XXXXXX'`, replacing XXXXXX with your assigned license key.
>
> Alternatively, the use of MODELLER can be manually disabled in the workflow steps.

Mamba and Snakemake are used to setup and manage the dependencies and workflows in this project. Conda must also be available for the inital setup. Run:

```
snakemake environment --cores 1
```

To create the mamba environment (will be located at ./.venv). Activate it by running:

```
mamba activate ./.venv/
```

### Creating AlphaFold containers

Two containers need to be created and placed in a "containers" directory (`mkdir containers`) at the project root to run the AlphaFold predictions in the workflow. These include:

* AlphaFold2 container: navigate from the project root into the TCRmodel2 directory (`cd third_party/tcrmodel2`) and build the container using `apptainer build alphafold2.sif singularity/tcrmodel2_singularity_cuda12-2.def`

* AlphaFold3 container: navigate into the AlphaFold 3 directory (`cd third_party/AlphaFold3`) and follow the [instructions here](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#running-using-singularity-instead-of-docker).

Once built, these containers need to be moved into the "containers" directory (`mv alphafold<#>.sif ../../containers`).

### Installing for Development

To install the required dependencies for development, run the above and additionally run the command:

```
pip install -e '.[develop]'
```

This will install the required development dependencies. To verfiy the installation, run:

snakemake test --cores 1

## Database Notes

Several sequence databases were manually downloaded to make up the datasets for this project:

1. IEDB
    * URL: https://www.iedb.org/
    * Date: 2025/01/14
    * Selection Criteria:
      * Assay: T cell, MHC Ligand, Outcome - Positive
      * Epitope: Any
      * MHC Restriction: Any
      * Host: Any
      * Disease: Any
      * Reference: Any
    * Export T Cell Receptors with Epitopes (single headers as csv)
    * Result: 224498 sequences downloaded

2. ITRAP
    * URL: https://doi.org/10.11583/DTU.22645342.v1
    * Date: 2025/01/14
    * Complete download of TCR (highest quality) dataset (2833 sequences)

3. McPas-TCR
    * URL: https://friedmanlab.weizmann.ac.il/McPAS-TCR/
    * Date: 2025/01/10
    * Download complete database (40779 sequences)

4. VDJdb
    * URL: https://vdjdb.cdr3.net/
    * Date: 2025/01/10
    * Selection Criteria:
      * CDR3 - Species: Human, Monkey, Mouse
      * CDR3 - Gene (chain): TRA, TRB
      * MHC - Class: MHCI, MHCII
      * Meta - Assay Type: Multimer sorting, Culture-based, Other
      * Meta - Sequencing: Sanger, High-throughput, Single-cell
      * Meta - Spurious CDR3: Include non-canonical, Include unmapped V/J
    * Result: All sequences downloaded (120849 of 120849)
