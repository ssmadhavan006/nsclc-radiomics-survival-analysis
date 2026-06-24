# Pipeline Reproducibility Guide

This document describes the environment specifications, package dependencies, CLI commands, and expected outputs required to reproduce the entire pipeline.

## 1. System Requirements and Runtime Info

*   **Operating System**: Windows / Linux / macOS (Windows 11 tested).
*   **Python Version**: `Python 3.9` (Recommended. PyRadiomics uses compiled C-extensions that may fail to compile or run on Python 3.10+ without visual studio build tools).
*   **Typical Runtime**:
    *   *Stage 2 (Ingestion)*: $\approx$ 15 seconds.
    *   *Stage 3 & 4 (Preprocessing & Feature Extraction)*: $\approx$ 20 minutes (across 20 CPU threads). Checks checkpoints and skips already processed patients.
    *   *Stage 5 & 6 (Cleaning & Statistics)*: $\approx$ 10 seconds.
    *   *Stage 7 (Survival)*: $\approx$ 45 seconds (including 1000-resample bootstrapping).
    *   *Total End-to-End Runtime*: $\approx$ 21 minutes (first run), under 2 minutes (subsequent runs due to checkpointing).

---

## 2. Recreating the Environment

Ensure `uv` is installed, then run:

```powershell
# 1. Create a virtual environment with Python 3.9
uv venv --python 3.9

# 2. Activate the virtual environment
# Windows (PowerShell):
.venv\Scripts\activate
# Linux/macOS:
source .venv/bin/activate

# 3. Install pinned dependencies
uv pip install -r requirements.txt
```

---

## 3. Pipeline Execution Commands

Execute the pipeline stages using the master CLI tool [run_pipeline.py](file:///d:/Coding/radiomics/run_pipeline.py):

```powershell
# Run the entire pipeline (Ingestion -> Extraction -> Analysis -> Survival)
python run_pipeline.py --stage all

# Run Ingestion Only (Parses dataset, creates manifest, logs failed_cases.csv)
python run_pipeline.py --stage ingestion

# Run Parallelized Feature Extraction Only (Skips already processed temp files)
python run_pipeline.py --stage extraction

# Run Data Cleaning & Statistics Only (Drops redundant features, creates PCA & Heatmap)
python run_pipeline.py --stage analysis

# Run Survival Modeling Only (Trains CPH, runs bootstrap validation, saves joblib models)
python run_pipeline.py --stage survival

# Launch the Gradio web dashboard demonstrator
python gradio_app.py
```

---

## 4. Expected Outputs and Artifacts

Upon successful execution of `python run_pipeline.py --stage all`, the following files will be generated in `outputs/`:

*   **Logs**:
    *   `outputs/logs/data_manifest.csv`: Status of all 422 patients scanned.
    *   `outputs/logs/failed_cases.csv`: Contains patient ID, pipeline stage, exclusion reason, and timestamp (e.g. `LUNG1-128` missing segmentation).
*   **Features**:
    *   `outputs/features/raw_features_all_patients.csv`: Combined feature matrix of 421 patients and 889 features.
    *   `outputs/features/cleaned_feature_matrix.csv`: Variance & Spearman correlation cleaned matrix of 421 patients and 190 features.
    *   `outputs/features/pyradiomics_params.yaml`: Exact settings utilized by the PyRadiomics extractor.
    *   `outputs/features/scaler.joblib`: Serialized `StandardScaler` fitted on the 19 selected features.
    *   `outputs/features/model_radiomics.joblib`: Serialized `CoxPHFitter` signature model.
    *   `outputs/features/model_combined.joblib`: Serialized `CoxPHFitter` combined clinical-radiomic model.
*   **Tables**:
    *   `outputs/tables/clinical_summary.csv`: Summary statistics of patient demographics.
    *   `outputs/tables/univariate_associations.csv`: Spearman/MW U/KW test p-values, FDR corrections, and effect sizes.
    *   `outputs/tables/cox_univariate_results.csv`: Standardized hazard ratios and C-index for each feature.
    *   `outputs/tables/cox_multivariate_results.csv`: Combined model coefficients, Wald CIs, bootstrap CIs, and p-values.
    *   `outputs/tables/model_comparison.csv`: Model comparison C-indices and 95% bootstrap intervals.
    *   `outputs/tables/feature_meanings.csv`: Interpretation mappings for the 19 selected features.
*   **Figures**:
    *   `outputs/figures/pca_plots.png`, `correlation_heatmap.png`, `feature_boxplots_survival.png`, `km_overall_survival.png`, `km_survival_by_stage.png`, `km_survival_by_histology.png`, `km_radiomic_signature.png`, `cox_multivariate_forest.png`, `time_dependent_roc.png`, `calibration_plot_3yr.png`.
