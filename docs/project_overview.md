# Project Overview: Quantitative Radiomic Characterization of NSCLC Tumors

## 1. Problem Statement
Lung cancer remains the leading cause of oncology-related mortality worldwide, with Non-Small Cell Lung Cancer (NSCLC) representing approximately 85% of all diagnoses. NSCLC is characterized by significant intratumoral spatial heterogeneity. Traditional clinical staging schemes—such as the tumor-node-metastasis (TNM) classification system—are coarse and fail to capture sub-visual microenvironmental variations. Consequently, patients with identical staging profiles often exhibit divergent therapeutic responses and survival outcomes.

## 2. Motivation
Clinical assessment of tumor morphology on Computed Tomography (CT) scans has historically relied on subjective qualitative visual interpretation (e.g., assessing borders, necrosis, and cavitation) or simple 1D/2D diameter measurements (RECIST criteria). 
Quantitative Radiomics overcomes these limitations by:
1.  **Extracting High-Dimensional Phenotypes**: Calculating hundreds of mathematical metrics (including volume, shape descriptors, first-order intensity distributions, and advanced wavelet textures).
2.  **Mitigating Inter-Observer Bias**: Providing objective, reproducible measurements computed programmatically.
3.  **Capturing Spatial Heterogeneity**: Assessing spatial patterns of voxel intensities, which reflect underlying biological factors like tumor vascularity, hypoxia, cellular density, and necrosis.

## 3. Dataset (TCIA NSCLC-Radiomics Lung1)
This research utilizes the public **NSCLC-Radiomics** (Lung1) dataset hosted on The Cancer Imaging Archive (TCIA). The cohort includes:
*   Pre-treatment chest CT scans of **422 patients** diagnosed with stage I-IIIb NSCLC.
*   Manual segmentations of the Gross Tumor Volume (GTV) outlined by radiation oncologists.
*   Complete clinical covariates (Age, Gender, Overall Stage, Histology, Overall Survival Time, and Vital Status).

## 4. Research Objectives
1.  **Standardized Preprocessing**: Implement coordinate-aligned 3D DICOM image and segmentation matching, intensity clipping, and isotropic resampling ($1.0 \times 1.0 \times 1.0\text{ mm}^3$) according to Image Biomarker Standardization Initiative (IBSI) guidelines.
2.  **Feature Extraction**: Extract 889 quantitative shape, first-order, and wavelet-filtered texture features from the GTV.
3.  **Univariate Association Checking**: Perform statistical association tests to link radiomics phenotypes with overall clinical stage and histological subtypes, applying Benjamini-Hochberg False Discovery Rate (FDR) corrections.
4.  **Prognostic Signature Modeling**: Build a robust, regularized LASSO-Cox Proportional Hazards model using cross-validation to construct a Prognostic Score signature.
5.  **Validation**: Validate model performance using 1000-resample bootstrapping, time-dependent Receiver Operating Characteristic (ROC) curves, and calibration analysis.

## 5. Summary of Contributions
*   Created a fully reproducible, modular pipeline that automates all stages from raw DICOM coordinates alignment to predictive modeling.
*   Identified 113 radiomic features showing significant correlation with clinical staging after FDR correction.
*   Developed a 19-feature radiomics signature that serves as an independent prognostic marker, outperforming standard clinical variables (C-index: 0.6396 vs. 0.5483).
*   Built and deployed an interactive Gradio app enabling clinical users to overlay GTV segmentation contours on axial slices and simulate prognosis changes using feature sliders.
