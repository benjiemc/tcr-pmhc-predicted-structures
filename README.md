# TCR:pMHC Predicted Structures

This project creates synthetic TCR:pMHC structures from sequence data using state of the art deep learning models.

## Installation and Setup

To install this project, first clone the repo:

```
git clone git@github.com:benjiemc/tcr-pmhc-predicted-structures.git
cd tcr-pmhc-predicted-structures/
```
<!--TODO add the databases needed-->

The following envrionment variables need to be set during the initialisation stage:

- ALPHAFOLD3_WEIGHTS_DIR: Path to the AF3 weights directory
- ALPHAFOLD_DATA_DIR: Path to the AlphaFold3 data directory
- COLABFOLD_DATA_DIR: Path to the ColabFold data directory
- TCRMODEL2_ALPHAFOLD_DATA_DIR: Path to the TCRModel2 data directory

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