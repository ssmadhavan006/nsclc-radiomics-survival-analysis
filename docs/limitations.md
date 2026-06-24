# Research Limitations

This document details the methodological and dataset-specific limitations of the current study.

## 1. Single-Center Scanner Protocol Bias
*   **Limitation**: The TCIA NSCLC-Radiomics Lung1 cohort was collected from a single center (MAASTRO Clinic, Netherlands).
*   **Impact**: While the scanning protocols are highly uniform (which benefits model training), it introduces scan acquisition bias.
*   **Risk**: The extracted radiomic features might overfit to the specific CT scanner settings (e.g. slice thickness, reconstruction kernel, and scanner manufacturer) used in the Lung1 cohort, potentially degrading performance when tested on external multi-center datasets.

## 2. Lack of External Validation
*   **Limitation**: The prognostic signature was developed and validated solely within the Lung1 dataset (using cross-validation and bootstrapping).
*   **Impact**: The performance of the signature has not yet been validated on an independent external validation dataset.
*   **Risk**: The reported concordance index ($C = 0.6396$) may represent an optimistic estimate. External validation is required to prove generalizability.

## 3. Retrospective Clinical Endpoints & Confounders
*   **Limitation**: The clinical outcomes data lacks information on post-imaging therapies and clinical management details.
*   **Impact**: Important clinical confounders—such as variable chemotherapy regimens, immunotherapy, targeted biological agents, radiation doses, and patient comorbidities—were not controlled for.
*   **Risk**: Differences in overall survival may be partially driven by heterogeneous treatment regimens rather than tumor phenotype alone, which could confound the association between radiomics and outcomes.

## 4. Tumor Segment Variability & Spatial Limits
*   **Limitation**: Segmentations were manually contoured by radiation oncologists.
*   **Impact**: Manual contouring is subject to inter-observer variability, which can introduce noise into shape and border-based radiomics metrics.
*   **Risk**: Minor variations in contour boundaries can alter texture features (especially high-frequency wavelet bands). Although coordinate validation was performed, segment outline variability remains an uncontrolled source of noise.
*   *Note*: The pipeline enforces a minimum threshold of 50 voxels, meaning small nodules are excluded from analysis.
