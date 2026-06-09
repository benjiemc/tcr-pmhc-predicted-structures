# TCR:pMHC Predicted Structures

This repository contains the code, scripts, and analysis to accompany the work "Characterising AlphaFold 3's ability to predict T cell antigen specificity". The analysis investigates of state-of-the-art deep learning models for structure prediction can be used to predict T cell antigen specificity.

## Overview

The results of the analysis are presented in the following notebooks:

| Notebook Name | Description | Associated Figures in Manuscript |
| ------------- | ----------- | -------------------------------- |
| [ablations_of_alphafold3.ipynb](notebooks/ablations_of_alphafold3.ipynb) | Ablation study of the effect of unpaired MSAs, paired MSAs, and templates on AlphaFold 3's performance | Figure 3 and S3 |
| [benchmark_alphafold3_against_sequence_models.ipynb](notebooks/benchmark_alphafold3_against_sequence_models.ipynb) | Comparison of AlphaFold 3's ability to predict T cell antigen specificity to NetTCR and TCRStructMap |  |
| [benchmark_structure_prediction_methods.ipynb](notebooks/benchmark_structure_prediction_methods.ipynb) | Benchmark of AlphaFold 2, AlphaFold 3, Chai-1, and Boltz-2 at predicting T cell antigen specificity and TCR:pMHC complex structures | Figure 2 and S2 |
| [c259_mutational_scanning.ipynb](notebooks/c259_mutational_scanning.ipynb) | Correlation anaylsis between the AlphaFold 3-derived binding scores and the affinity data of the c259 TCR binding to a peptide mutational scanning library | Figure S8 |
| [coeliac_antigen_specificity_with_alphafold3.ipynb](notebooks/coeliac_antigen_specificity_with_alphafold3.ipynb) | Assessing the ability of AlphaFold 3 to predict the specificity of coeliac TCRs binding to gluten peptides | |
| [coeliac_mutagenesis_analysis.ipynb](notebooks/coeliac_mutagenesis_analysis.ipynb) | Analysis of the effects of mutations to the 572G TCR binding score compared to affinity data | Figure 7 |
| [compare_msa_generation.ipynb](notebooks/compare_msa_generation.ipynb) | Comparison of the different MSA strategies run time, ability to predict T cell antigen specificity, and ability to predict TCR:pMHC complex structures | Figure 1 and S1 |
| [compare_structural_features_to_confidence.ipynb](notebook/compare_structural_features_to_confidence.ipynb) | Correlation analysis of the AlphaFold 3 confidence scores to structural features of TCR:pMHC complexes | Figure 4 and S4, Table 1|
| [dash_mutagenesis_analysis.ipynb](notebooks/dash_mutagenesis_analysis.ipynb) | Analysis of the effect of mutations on AlphaFold confidence scores for TCRs reported in Dash et al. binding the BMLF epitope | |
| [evaluating_antigen_specificity_model.ipynb](notebooks/evaluating_antigen_specificity_model.ipynb) | Analysing the performace of PAE Aggregator Scores at predicting T cell antigen specificity | Figures S5 |
| [exploring_coeliac_tcrs.ipynb](notebooks/exploring_coeliac_tcrs.ipynb) | Data exploration and clustering of TCRs reacting to coeliac antigens | |
| [exploring_dash_dataset.ipynb](notebooks/exploring_dash_dataset.ipynb) | Data exploration and clustering of TCRs reacting to antigens reported by Dash et al. | Figure 6, S6, and S7 |
| [exploring_specificity_prediction_normalisation.ipynb](notebooks/exploring_specificity_prediction_normalisation.ipynb) | Looking at the impact of the PAE Aggregator Scores on normalising TCR:pMHC specificity predictions | Figure 5 |
| [investigate_specificity_information_in_af3_scores.ipynb](notebooks/investigate_specificity_information_in_af3_scores.ipynb) | Comparison of different strategies to predict T cell antigen specificity using AlphaFold 3 confidence scores | |
| [visualise_sequence_data.ipynb](notebooks/visualise_sequence_data.ipynb) | Visualisations of the sequence data used in this project. | |
| [visualise_structure_data.ipynb](notebooks/visualise_structure_data.ipynb) | Visualisations of the structure data used in this project. | |

To run the analysis, follow the instructions below and use the [Snakemake workflow](workflow/Snakefile).

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

> ### Using NetMHCpan
> A small part of the analysis uses NetMHCpan version 4.2 to screen out non-presentable peptides.
> In order to run this analysis, NetMHCpan 4.2 needs to be downloaded and installed in the environment following the instructions [here](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/).

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

```
snakemake test --cores 1
```

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
