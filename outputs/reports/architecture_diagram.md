# Quantitative Radiomics Pipeline Architecture Diagram

This document illustrates the end-to-end quantitative radiomics pipeline for characterization of non-small cell lung cancer (NSCLC) tumors and overall survival outcomes prediction.

## Research-Grade Architecture Diagram

Below is the generated high-impact pipeline diagram:

![End-to-End Quantitative Radiomics Pipeline Diagram](file:///C:/Users/Admin/.gemini/antigravity/brain/3ae17894-4049-43e6-8f28-e1bae0b90fff/architecture_diagram_1782315322768.png)

---

## Detailed Pipeline Stages

### 1. Data Ingestion & Quality Control
* **Inputs**: Clinical metadata spreadsheet (CSV) and raw 3D DICOM CT scans accompanied by GTV primary tumor RT-struct segmentations.
* **Checks**: Discovers missing volumes/masks (e.g. `LUNG1-128`) and logs exclusions to `failed_cases.csv`. Sorts DICOM slices by Z-coordinate.

### 2. Image Preprocessing
* **Clipping**: CT voxel intensities clipped to $[-1000, 400]$ HU to isolate soft-tissue and lung densities while discarding noise.
* **Resampling**: Isotropic voxel spacing standardization to $1.0 \times 1.0 \times 1.0\text{ mm}^3$ using standard interpolation (B-Spline for CT volumes, Nearest Neighbor for segmentation masks).
* **Overlap QC**: Aborts feature extraction if CT dimensions, origins, or directions mismatch the segmentation mask coordinate grid.

### 3. High-Throughput Radiomic Feature Extraction
* **Core Extractor**: Runs parallel PyRadiomics engine using multi-threaded Loky backend (20 CPU threads) with temporary checkpointing.
* **Feature Vector**: Extracts 889 features per patient:
  * 14 3D Shape features
  * 18 First-Order statistics features
  * 75 Texture features (GLCM, GLSZM, GLRLM, NGTDM, GLDM)
  * 782 Wavelet decompositions (decomposed across LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH sub-bands)

### 4. Dimensionality Reduction & Cleaning
* **Variance Threshold**: Filters out 289 low-variance features ($Var < 0.01$) that carry minimal spatial texture information.
* **Spearman Redundancy**: Drops 372 collinear features with pairwise Spearman correlation coefficients $|\rho| > 0.95$, retaining 190 independent phenotypes.

### 5. Statistical Association Testing
* **Staging Association**: Mann-Whitney U, Spearman $\rho$, and Kruskal-Wallis tests corrected via Benjamini-Hochberg False Discovery Rate (FDR) at $\alpha = 0.05$. Selects 113 stage-associated phenotypes.

### 6. Signature Modeling & Validation
* **Feature Selection**: 5-fold cross-validated LASSO-penalized Cox Proportional Hazards regression shrinks coefficients to isolate a 19-feature radiomic signature.
* **Prognostic Score**: Signature partial hazard is scaled against cohort median.
* **Bootstrap Validation**: 1000 bootstrap resamples are generated to calculate 95% Confidence Intervals for C-index and Wald multivariate Hazard Ratios.
* **Performance Evaluation**: Generates time-dependent ROC curve AUC values (1, 3, 5-year overall survival) and calibration quintile alignments.
