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

To create the mamba environment (will be located at ./.venv).

### Installing for Development

To install the required dependencies for development, run the above and additionally run the command:

```
pip install -e '.[develop]'
```

This will install the required development dependencies. To verfiy the installation, run:

snakemake test --cores 1