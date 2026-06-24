# Data Dictionary: Clinical Variables and Feature Groups

This document defines the schema of the clinical outcomes dataset and the categorization of the quantitative radiomic features.

## 1. Clinical Outcome Variables

The clinical database ([NSCLC-Radiomics-Lung1.clinical-version3-Oct-2019.csv](file:///d:/Coding/radiomics/dataset/NSCLC-Radiomics-Lung1.clinical-version3-Oct-2019.csv)) contains the following variables:

| Variable Name | Data Type | Units / Range | Description |
| :--- | :--- | :--- | :--- |
| **PatientID** | String | `LUNG1-001` to `LUNG1-422` | Unique patient identification code. |
| **age** | Float | Years (Range: 35.0 - 89.0) | Age of the patient at the time of diagnosis. |
| **clinical.T.Stage** | String | `T1`, `T2`, `T3`, `T4`, `Tx` | Clinical classification of the primary tumor size and extent. |
| **Clinical.N.Stage** | String | `N0`, `N1`, `N2`, `N3`, `Nx` | Clinical classification of regional lymph node involvement. |
| **Clinical.M.Stage** | String | `M0`, `M1`, `Mx` | Clinical classification of distant metastasis. |
| **Overall.Stage** | String | `I`, `II`, `IIIa`, `IIIb`, `Unknown` | Overall clinical stage grouping combining TNM status. |
| **Histology** | String | Squamous Cell Ca, Adeno Ca, Large Cell, etc. | Histological classification of the tumor cell type. |
| **gender** | String | `male`, `female` | Biological gender of the patient. |
| **Survival.time** | Float | Days (Range: 5.0 - 3683.0) | Overall survival time measured from diagnosis to death or last follow-up. |
| **deadstatus.event** | Integer | `0` (Censored), `1` (Deceased) | Vital status endpoint indicator. |

---

## 2. Radiomic Feature Categories

The raw features matrix contains 889 features, grouped into the following mathematical categories:

| Feature Class | Raw Count | Wavelet Count | Description |
| :--- | :--- | :--- | :--- |
| **Shape descriptors** | 14 | 0 | Geometric descriptors representing volume, surface area, compactness, and elongation of the 3D GTV mask. Shape features are not calculated on wavelet-filtered volumes. |
| **First-Order statistics** | 18 | 144 | Voxel intensity statistics (e.g. mean, median, skewness, kurtosis, range, energy) describing the distribution inside the GTV. |
| **GLCM** | 24 | 192 | Gray Level Co-occurrence Matrix features representing spatial relationships between pairs of voxels with specific gray-level values. |
| **GLRLM** | 16 | 128 | Gray Level Run Length Matrix features describing the distribution of consecutive voxels with the same gray-level value in a given direction. |
| **GLSZM** | 16 | 128 | Gray Level Size Zone Matrix features describing the size and intensity distribution of connected 3D gray-level zones. |
| **GLDM** | 14 | 112 | Gray Level Dependence Matrix features representing dependencies of voxel intensities on neighbors within a defined distance. |
| **NGTDM** | 5 | 40 | Neighborhood Gray Tone Difference Matrix features capturing local texture properties like coarseness, contrast, business, complexity, and strength. |
| **Total Features** | **145** | **744** | **Sum Total: 889 features extracted per patient scan.** |
