# Feature Generation Pipeline

This directory contains the MATLAB and Python scripts used to generate the feature tables used in the machine learning analyses presented in the accompanying study.

Feature generation was performed separately for the original *Saccharomyces cerevisiae* genome and for the randomized genomes used in the validation experiments.

---

## Common dataset preparation

The following scripts generate the sequence datasets used by both the original and randomized genome pipelines.

- `gen_Y_genes.m` – generates the initial gene dataset from the reference genome and annotation files.
- `orig_3utrs.m` – reconstructs the original 3′UTR sequences and generates `Y_genes_orig3utrs.mat`, which is used as the input for downstream feature generation.

---

## Original genome pipeline

The original genome pipeline generates the novel features introduced in this study.
Previously published feature values were obtained from the published dataset described in the study and therefore are not recalculated for the original genome.

### 1. Weight calculation

The following scripts generate the weight tables required for CAI and tAI calculation.

- `CAI_weights.m`
- `tAI_weights.m`

### 2. Novel feature generation

The following scripts calculate the novel codon-related features introduced in this study.

- `CAI_feature.m`
- `tAI_feature.m`
- `codons_features.m`
- `tRNA_feature.m`

### 3. Stop-window MFE calculation

The stop-window MFE feature is generated in three steps.

1. Generate FASTA files
   - `fastas_stop_win.m`

2. Calculate minimum free energy using the ViennaRNA package
   - `main.py`

3. Parse the output into a MATLAB feature table
   - `stop_win_mfe.m`

### 4. Final feature table

- `reorg_features_table.m` aligns the generated feature tables.
- `gen_XYtbls_wt25.m` combines the novel features with the previously published feature set and exports the final feature matrix used for machine learning analyses.

---

## Randomized genome pipeline

The randomized genome pipeline generates randomized coding sequences and recalculates **both** the previously published features and the novel features introduced in this study.

### 1. Randomized genome generation

- `gen_Y_rand.m`
- `perm_nt.m`
- `perm_codons.m`

### 2. Previously published feature generation

The following scripts calculate the previously published features for the randomized genomes.

- `Y_features.m`
- `utr3stop_feature.m`

The original 3′UTR MFE feature is generated in three steps.

1. Generate FASTA files
   - `orig3utrs_fastas.m`

2. Calculate minimum free energy using the ViennaRNA package
   - `main.py`

3. Parse the output into a MATLAB feature table
   - `orig3utrs_mfe_feature.m`

### 3. Novel feature generation

The following scripts calculate the novel codon-related features.

- `CAI_feature.m`
- `tAI_feature.m`
- `codons_features.m`
- `tRNA_feature.m`

The stop-window MFE feature is generated using the same workflow as for the original genome:

1. `fastas_stop_win.m`
2. `main.py`
3. `stop_win_mfe.m`

### 4. Final feature tables

- `reorg_features_table.m` aligns the generated feature tables.
- `gen_Xtbls.m` exports the randomized feature matrices used for the validation analyses.