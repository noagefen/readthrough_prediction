# Determinants of Translation Readthrough in Eukaryotes

This repository contains the feature-generation and machine-learning code for the accompanying study of translation readthrough in *Saccharomyces cerevisiae*.

## Current analysis

The revised analysis is available in two equivalent forms:

- `notebooks/compare-models-nested-cv.ipynb` is the clean, upload-ready notebook.
- `analysis/nested_cv_analysis.py` contains the same executable analysis as a Python script.

Both versions compare the previously tested feature set with the full feature set introduced in this study. Numeric features are standardized within each training fold, categorical features are one-hot encoded, and hyperparameters are selected by inner cross-validation without using the held-out evaluation data. The randomized-transcript analysis uses matched out-of-fold predictions for each native gene and its 20 randomized counterparts.

The pipeline generates the revised model-performance summaries, Figures 3, 4, and 6, Supplementary Figures S1-S3, and the data used for Supplementary Tables S3-S4. Generated files are written to `nested_cv_results/` and are excluded from version control.

## Running the analysis

Install the exact package versions:

```bash
python -m pip install -r requirements-nested-cv.txt
```

Place the model inputs in the following directories at the repository root:

```text
WT25/
  Xtbl_wt25.csv
  Ytbl_wt25.csv
  train_inds_tbl.csv
  test_inds_tbl.csv
data/
  Xtbl_rand_1.csv
  ...
  Xtbl_rand_20.csv
```

Then run either the notebook or:

```bash
python analysis/nested_cv_analysis.py
```

The code also recognizes the Kaggle datasets [`wt25-dataset`](https://www.kaggle.com/datasets/noagef/wt25-dataset) and [`randgenome-x`](https://www.kaggle.com/datasets/noagef/randgenome-x), which contain the exact input tables used for the reported analysis. The input tables are not duplicated in this repository because of their size.

## Repository structure

- `analysis/` contains the current reproducible analysis script.
- `notebooks/compare-models-nested-cv.ipynb` contains the current notebook.
- `feature_extraction/orig_genome/` contains MATLAB and Python scripts for the native-transcript feature table.
- `feature_extraction/rand_genomes/` contains scripts for the randomized-transcript feature tables.
- `notebooks/compare-models.ipynb` and `notebooks/randomized-genome.ipynb` are retained as the originally submitted analysis notebooks.

The source datasets and feature definitions are described in the manuscript and its data-availability statement.
