# Research Methodology: Quantitative Radiomics Pipeline

This document describes the methodological framework, preprocessing algorithms, feature selection steps, and survival modeling techniques implemented in the pipeline.

## 1. Patient Cohort & Ingestion QC
*   **Initial Scanned Cohort**: 422 patients from the TCIA NSCLC-Radiomics Lung1 dataset.
*   **Quality Control (QC)**: Validated for the presence of a complete CT volume, a valid coordinate-overlapping segmentation mask, and complete outcome data.
*   **Exclusion**: Patient `LUNG1-128` was excluded due to a missing segmentation series (RTSTRUCT/SEG candidates absent) and logged in `failed_cases.csv`.
*   **Final Study Cohort**: **421 patients** successfully processed.

## 2. Preprocessing & Image Standardization
To ensure feature reproducibility and satisfy Image Biomarker Standardization Initiative (IBSI) guidelines:
1.  **Slice Ordering**: Slices within each CT scan were explicitly sorted ascending by `ImagePositionPatient[2]` (Z-coordinate) to resolve variable slice order anomalies.
2.  **Coordinates Alignment**: Computed the 3D coordinate transformation matrices for the multi-frame DICOM SEG mask to match the primary CT's spatial grid:
    *   Target segment name: "Neoplasm, Primary", "Primary", or "Neoplasm".
    *   2D frame pixels mapped to CT Z-coordinates and stacked into a 3D binary volume matching CT spatial information.
3.  **Intensity Normalization (HU Scaling)**: Clipped CT voxel intensities to Hounsfield Unit (HU) range $[-1000, 400]$ to remove background air ($<-1000$) and high-density structures like bone or metal artifacts ($>400$).
4.  **Isotropic Resampling**: Resampled CT volumes and binary segmentations to a standard $1.0 \times 1.0 \times 1.0\text{ mm}^3$ voxel spacing using SimpleITK:
    *   CT image: `sitkBSpline` interpolation.
    *   Binary mask: `sitkNearestNeighbor` interpolation (preserving binary labels, followed by a threshold filter at 1.0).
5.  **Spatial Alignment QC Verification**: Spacings, origins, direction cosines, and grid dimensions were verified. The GTV mask was validated to contain $\ge 50$ voxels.

## 3. High-Throughput Feature Extraction
Using PyRadiomics, 889 features were extracted from the resampled GTV (within the clipped CT volume):
*   **Shape-based Features** (14)
*   **First-Order Statistics** (18)
*   **Texture Features** (75): GLCM (24), GLRLM (16), GLSZM (16), GLDM (14), and NGTDM (5).
*   **Wavelet Features** (744): 8 decomposition levels (HHH, HHL, HLH, HLL, LHH, LHL, LLH, LLL).
*   *Parameters*: Fixed bin discretization width of 25 HU, normalization scale of 100.

## 4. Feature Cleaning & Engineering
To mitigate multi-collinearity and overfitting, the feature space was reduced from 889 to 190 dimensions:
1.  **Variance Filtering**: Dropped 289 features with variance $< 0.01$, eliminating near-constant features.
2.  **Correlation Filtering**: Calculated pairwise Spearman rank correlation coefficient ($\rho$) matrix of the remaining features. For features with $|\rho| > 0.95$, the average correlation to all other features was calculated. The feature with the higher average correlation was dropped (372 features dropped).

## 5. Prognostic Signature Construction
1.  **Data Partitioning**: Standard Scaling was fit on training folds and applied to test folds to prevent data leakage (Rule 28).
2.  **LASSO-Cox Model Tuning**: Nested 5-fold cross-validation was used to tune the L1 penalizer parameter ($p \in [0.01, 0.05, 0.1, 0.2, 0.3]$).
3.  **Signature Extraction**: At the optimal regularizer penalizer ($0.05$), the LASSO penalty forced non-predictive feature coefficients to zero, selecting **19 features** to construct the Prognostic Score:
    \[\text{Prognostic Score} = \exp\left(\sum_{i=1}^{19} \beta_i \cdot Z_i\right)\]
    where $Z_i$ is the standardized feature value and $\beta_i$ is the corresponding Cox regression coefficient.

## 6. Validation Framework
*   **Model Comparison**: Evaluated and compared three Cox models: Clinical variables only, Radiomics Prognostic Score only, and Combined.
*   **Bootstrap Validation**: 1000-resample bootstrapping was performed on the cohort to construct empirical 95% Confidence Intervals for C-index values and multivariate Hazard Ratios.
*   **Time-Dependent Discrimination**: Calculated Area Under the Curve (AUC) for predicting overall survival at 1-year, 3-year, and 5-year timepoints.
*   **Calibration Testing**: Binned predicted survival probabilities into quintiles at 3-years and plotted them against observed Kaplan-Meier survival curves.
