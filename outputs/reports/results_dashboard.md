# Quantitative Radiomics Validation Dashboard (TCIA NSCLC Lung1)

This dashboard consolidates the verified quantitative metrics and statistical results of the end-to-end NSCLC-Radiomics pipeline (421 patients).

---

## 1. Key Performance Metrics

| Metric | Value | 95% Confidence Interval | Source File |
|---|---|---|---|
| **Clinical-Only C-Index (Model A)** | `0.5483` | `[0.5182 - 0.5854]` (Bootstrap) | `outputs/tables/model_comparison.csv` |
| **Radiomics-Only C-Index (Model B)** | `0.6396` | `[0.6084 - 0.6693]` (Bootstrap) | `outputs/tables/model_comparison.csv` |
| **Combined Model C-Index (Model C)** | `0.6328` | `[0.6079 - 0.6624]` (Bootstrap) | `outputs/tables/model_comparison.csv` |
| **1-Year Survival AUC** | `0.700` | N/A | Dynamic (Time-Dependent ROC) |
| **3-Year Survival AUC** | `0.730` | N/A | Dynamic (Time-Dependent ROC) |
| **5-Year Survival AUC** | `0.714` | N/A | Dynamic (Time-Dependent ROC) |
| **Prognostic Score Hazard Ratio** | `2.554` | `[2.194 - 3.811]` (Bootstrap) | `outputs/tables/cox_multivariate_results.csv` |
| **Age Hazard Ratio** | `1.014` | `[1.002 - 1.028]` (Bootstrap) | `outputs/tables/cox_multivariate_results.csv` |
| **Stage (III vs I/II) Hazard Ratio** | `0.944` | `[0.767 - 1.180]` (Bootstrap) | `outputs/tables/cox_multivariate_results.csv` |
| **Gender (Male vs Female) HR** | `1.175` | `[0.927 - 1.499]` (Bootstrap) | `outputs/tables/cox_multivariate_results.csv` |

---

## 2. Ingestion & Quality Control Summary

* **Scanned Patients**: 422
* **Valid Patients**: 421
* **Excluded Patients**: 1 (Patient `LUNG1-128` - missing GTV segmentation; logged in `outputs/logs/failed_cases.csv`)
* **HU Clipping Bounds**: `[-1000, 400]` Hounsfield Units
* **Target Resampling Voxel Size**: `1.0 x 1.0 x 1.0 mm^3` (Isotropic, sitkBSpline/sitkNearestNeighbor)

---

## 3. Feature Pipeline Dimensions

* **Extracted Raw Features**: 889
* **Low Variance Removed (Var < 0.01)**: 289 features
* **Highly Redundant Removed (|rho| > 0.95)**: 372 features
* **Cleaned Features Matrix**: 190 features (`outputs/features/cleaned_feature_matrix.csv`)
* **LASSO-Cox Selected Features**: 19 features (`outputs/tables/feature_meanings.csv`)

---

## 4. Statistical Association Highlights

* **Stage Association**: **113 / 190 features** are statistically significant under Spearman correlation with Benjamini-Hochberg FDR correction ($\alpha = 0.05$).
* **Histology Association**: No individual features showed statistically significant association after multiple-testing correction.
* **Survival Status Association**: No individual features showed statistically significant univariate association after multiple-testing correction, indicating the necessity of multivariate signature modeling (LASSO-Cox).

---

## 5. Cohort Clinical Characteristics

### Overall Stage Distribution ($N=421$)
- **Stage IIIb**: 176 (41.8%)
- **Stage IIIa**: 111 (26.4%)
- **Stage I**: 93 (22.1%)
- **Stage II**: 40 (9.5%)
- **Stage Unknown**: 1 (0.2%)

### Histological Subgroup Distribution ($N=421$)
- **Squamous Cell Carcinoma**: 152 (36.1%)
- **Large Cell**: 114 (27.1%)
- **Nos (Not Otherwise Specified)**: 62 (14.7%)
- **Adenocarcinoma**: 51 (12.1%)
- **Unknown/Missing**: 42 (10.0%)
