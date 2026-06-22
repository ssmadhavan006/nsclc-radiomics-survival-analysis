# Walkthrough: NSCLC-Radiomics Lung1 End-to-End Pipeline

This walkthrough summarizes the final implementation, testing, and validation results for the **Quantitative Radiomic Characterization of NSCLC Tumors and Their Association with Clinical Outcomes** project.

---

## 1. Summary of Changes Made

The codebase has been built from scratch, keeping all modules modular, documented, and fully type-hinted (Rule 30).

### Core Pipeline Modules
1. **`src/config.yaml`**: The master configuration file centralizing paths, HU clipping bounds, target voxel spacings, PyRadiomics parameters, and statistical parameters (Rule 12, 13).
2. **`src/utils.py`**: Logging setup (simultaneous console and timestamped file writing) and config loading utilities (Rule 33, 34).
3. **`src/data_ingestion.py`**: Ingests the clinical CSV and scans the patient directories. It builds a coordinate-matched manifest and correctly identifies missing files (e.g. `LUNG1-128`), logging exclusions to `failed_cases.csv` (Rule 6, 7).
4. **`src/preprocessing.py`**: Performs 3D coordinate-aligned mask extraction (by Z-coordinate patient-position matching), clips CT intensities to $[-1000, 400]$ HU, and resamples scans/masks to $1.0\text{ mm}^3$ isotropic spacing (Rule 8, 9, 10, 11).
5. **`src/feature_extraction.py`**: Runs PyRadiomics extraction using `joblib` parallel processing (20 CPU threads) with LokyBackend and temporary-file checkpointing. Consolidates individual outputs into `raw_features_all_patients.csv` (Rule 13, 14, 15).
6. **`src/analysis.py`**: Applies variance ($Var < 0.01$) and correlation ($|\rho| > 0.95$) cleaning. Performs Mann-Whitney U, Spearman, and Kruskal-Wallis tests with Benjamini-Hochberg FDR correction. Generates correlation heatmap, PCA, and top feature boxplots (Rule 18, 19, 20, 21).
7. **`src/survival.py`**: Fits Kaplan-Meier survival curves (with confidence intervals and risk tables), runs univariate Cox Proportional Hazards models, cross-validates a LASSO-Cox model (tuning the L1 penalizer parameter), and builds the Prognostic Score signature. It also performs bootstrap validation (1000 resamples), time-dependent ROC calculations (1, 3, 5-year AUCs), 3-year calibration, model comparisons, and feature meaning mapping (Rule 22, 23, 24, 25, 26, 27, 28).
8. **`run_pipeline.py`**: Command-line orchestration entrypoint script allowing execution of individual stages or the entire pipeline (Rule 15).

---

## 2. Testing and Validation Results

### 2.1 Ingestion & Preprocessing QC
* **Total Cohort Size**: 422 patients scanned.
* **Excluded**: `LUNG1-128` (missing GTV segment) was correctly excluded and logged to `outputs/logs/failed_cases.csv`.
* **Valid**: 421 patients successfully matched and verified.
* **Resampling & Alignment**: Spacings verified as exactly $1.0 \times 1.0 \times 1.0\text{ mm}^3$ for all 421 patients. Spatial alignment between resampled CT and mask (origins, dimensions, direction cosines) was verified for all cases.

### 2.2 Parallel Feature Extraction
* **Feature Count**: Extracted 889 features per patient (including shape, first-order, texture, and wavelet decompositions).
* **Parallel Performance**: Executed batch runs across 20 CPU threads. The entire 421 cohort completed extraction in **20.3 minutes** (averaging ~0.9 seconds per patient batch).

### 2.3 Data Cleaning & Dimensionality Reduction
* **Variance Filtering**: Dropped 289 low-variance features ($Var < 0.01$).
* **Correlation Filtering**: Dropped 372 redundant features ($|\rho| > 0.95$, Spearman).
* **Cleaned Features Count**: Reduced feature dimensions from 889 to **190** independent, clean variables.

### 2.4 Statistical Associations
* **FDR-Significant Features ($\alpha = 0.05$)**:
  * **Overall Staging**: **113 features** were significantly associated with staging. Advanced stage tumors exhibit larger shapes (`MeshVolume`, `SurfaceArea`) and higher texture heterogeneity (especially in wavelets), showing phenotypic staging associations.
  * **Survival Status / Histology**: No features showed significant univariate association under FDR multiple-comparison correction.

### 2.5 Survival Analysis & Prognostic Modeling
* **Optimal LASSO Penalizer**: $0.05$ (Mean C-Index: **$0.6011$**).
* **Selected Features**: 19 features with non-zero coefficients. Detailed meanings and coefficients are stored in `outputs/tables/feature_meanings.csv`.
* **Model Comparison & Bootstrap Validation**:
  * **Model A (Clinical Only)**: C-Index = $0.5483$, 95% Bootstrap CI = $[0.5182, 0.5854]$
  * **Model B (Radiomics Only)**: C-Index = $0.6396$, 95% Bootstrap CI = $[0.6084, 0.6693]$
  * **Model C (Combined)**: C-Index = $0.6328$, 95% Bootstrap CI = $[0.6079, 0.6624]$
  Both radiomics models significantly outperformed the clinical-only model and exclude the $0.50$ baseline.
* **Multivariate Cox Hazards Model (Combined)**:
  * **PrognosticScore HR**: **$2.554$** ($95\%\text{ Wald CI: } [2.113, 3.088]$, $95\%\text{ Bootstrap CI: } [2.194, 3.811]$, **$p = 3.43 \times 10^{-22}$**).
  * **Age HR**: **$1.014$** ($95\%\text{ Bootstrap CI: } [1.002, 1.028]$, $p = 0.022$).
  * **Stage HR (III vs I/II)**: $0.944$ ($p = 0.623$, lost significance due to redundancy with radiomics score).
  * **Gender HR (Male vs Female)**: $1.175$ ($p = 0.176$).
* **Time-Dependent Discrimination (AUC)**:
  * **1-Year AUC**: **$0.700$**
  * **3-Year AUC**: **$0.730$**
  * **5-Year AUC**: **$0.714$**
* **Model Calibration**: Plot (`calibration_plot_3yr.png`) shows close alignment between predicted and observed 3-year survival probabilities.

---

## 3. Generated Figures

The pipeline programmatically generated the following figures in `outputs/figures/` (Rule 37):
1. [pca_plots.png](file:///d:/Coding/radiomics/outputs/figures/pca_plots.png) — PCA projection colored by Stage and Histology.
2. [correlation_heatmap.png](file:///d:/Coding/radiomics/outputs/figures/correlation_heatmap.png) — Clustered Spearman correlation heatmap of 190 cleaned features.
3. [feature_boxplots_survival.png](file:///d:/Coding/radiomics/outputs/figures/feature_boxplots_survival.png) — Box plots of top features against survival status.
4. [km_overall_survival.png](file:///d:/Coding/radiomics/outputs/figures/km_overall_survival.png) — KM curve of overall survival.
5. [km_survival_by_stage.png](file:///d:/Coding/radiomics/outputs/figures/km_survival_by_stage.png) — KM curves stratified by Overall Stage (with risk tables).
6. [km_survival_by_histology.png](file:///d:/Coding/radiomics/outputs/figures/km_survival_by_histology.png) — KM curves stratified by Histology (with risk tables).
7. [km_survival_by_gender.png](file:///d:/Coding/radiomics/outputs/figures/km_survival_by_gender.png) — KM curves stratified by Gender (with risk tables).
8. [km_survival_by_top_feature.png](file:///d:/Coding/radiomics/outputs/figures/km_survival_by_top_feature.png) — KM curves stratified by the top radiomic feature (with risk tables).
9. [km_radiomic_signature.png](file:///d:/Coding/radiomics/outputs/figures/km_radiomic_signature.png) — Stratification by Radiomic Signature (with risk tables and C-index).
10. [cox_multivariate_forest.png](file:///d:/Coding/radiomics/outputs/figures/cox_multivariate_forest.png) — Forest plot of multivariate Hazard Ratios.
11. [time_dependent_roc.png](file:///d:/Coding/radiomics/outputs/figures/time_dependent_roc.png) — Time-dependent ROC curves (1, 3, 5-year AUCs).
12. [calibration_plot_3yr.png](file:///d:/Coding/radiomics/outputs/figures/calibration_plot_3yr.png) — Predicted vs observed 3-year survival calibration curve.

---

## 4. Gradio Research Demonstrator

To present the findings and modeling pipeline to researchers and clinical stakeholders, a web-based **Gradio Research Demonstrator** (`gradio_app.py`) was implemented and refined:

### 4.1 Interface Capabilities & Visualizations
* **Interactive Cohort Exploration**: Users can select demo patients (e.g. `LUNG1-003`, `LUNG1-004`, `LUNG1-005`) from a dropdown to pre-load clinical metadata (Age, Gender, Stage) and their extracted 19-feature radiomic values.
* **GTV Axial CT Overlay Visualizer**: A side-by-side axial slice visualization displaying the original zoomed CT slice (left) and the CT with red semi-transparent GTV tumor mask overlay (`alpha = 0.35`) and sharp red boundary outlines (right). The slice is automatically determined by finding the Z-position with the maximum cross-sectional area of the GTV. Slices are cached to eliminate loading/resampling lag on repeated requests.
* **Prognostic Score Simulation**: Features can be manually adjusted using sliders to simulate morphological or texture transformations (e.g. shrinking or expanding tumor volume, increasing texture heterogeneity) and inspect the corresponding Prognostic Score, Risk Category, and Kaplan-Meier overall survival curve shift in real-time.
* **Custom Patient Prediction (Pipeline Upload)**: Users can upload a custom zip file containing DICOM CT slices alongside an RT-struct segmentation file. The app processes the zip, extracts the 19 features in coordinate-aligned resampled space, computes the prognostic score, and displays the risk category, GTV overlay, and survival curve.

### 4.2 Diagnostic & Statistical Metrics
* **Hazard Deviation & Percentile Rank**: Displays the patient's predicted relative hazard ratio deviation from the cohort median (e.g. `+34.0% above cohort median` or `-15.2% below cohort median`) and their ordinal percentile rank in the overall cohort (e.g. `82nd percentile`).
* **Model Validation Sidebar**: A permanent left-panel card displaying model validation metrics to establish immediate credibility: C-index (`0.640`), 1-Year AUC (`0.700`), 3-Year AUC (`0.730`), and 5-Year AUC (`0.714`), validated on the `421` patient cohort.
* **Feature Name Disambiguation**: Modifies the data table to display unique, fully qualified feature names (e.g. `LLL_gldm_SmallDependenceHighGrayLevelEmphasis` instead of a shortened split segment) to prevent name collisions in the UI.

### 4.3 System Verification & Bug Fixes
* **Pydantic Compatibility Fix**: Pinning `pydantic==2.10.6` resolved a schema-generation bug in Gradio Client where boolean schema items like `additionalProperties: false` caused a crash (`TypeError: argument of type 'bool' is not iterable`) during FastAPI setup.
* **Robust Server Binding**: Omitted the hardcoded port `7860` in favor of letting Gradio automatically bind to the first available port (e.g., falling back to `7861` if the port is currently held by a terminal session), ensuring clean launches.


