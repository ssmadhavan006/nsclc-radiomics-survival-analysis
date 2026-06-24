# Future Work Directions

This document outlines next-step research directions to expand, validate, and clinically translate the radiomics pipeline.

## 1. External Cohort Validation
*   **Objective**: Validate the 19-feature prognostic signature on external multi-center cohorts.
*   **Target Datasets**: Test on public datasets such as the **NSCLC-Radiomics-Interobserver** cohort or datasets from the Clinical Proteomic Tumor Analysis Consortium (CPTAC-LUAD/LSCC) to evaluate the signature's robustness to varying scanner protocols and reconstruction kernels.

## 2. Radiogenomics Integration
*   **Objective**: Integrate radiomics phenotypes with genomic, transcriptomic, and mutational datasets (e.g. EGFR, ALK, KRAS status).
*   **Implementation**: Link extracted GTV features with matching molecular profiles (e.g., using the TCIA/TCGA NSCLC dataset) to build multi-modal models that map radiographic texture to tumor biology.

## 3. Deep Survival Modeling
*   **Objective**: Evaluate deep-learning-based survival architectures to compare with classical Cox Proportional Hazards modeling.
*   **Implementation**: Train neural network architectures like **DeepSurv** or **Cox-nnet** on the raw 3D CT volumes and GTV masks to compare features learned by deep networks with classical handcrafted PyRadiomics features.

## 4. Explainable Radiomics (XAI)
*   **Objective**: Bridge the gap between mathematical radiomic metrics and physical histology.
*   **Implementation**: Use spatial visualization techniques (such as heatmaps mapping the intensity of specific wavelet features onto GTV coordinates) to allow clinicians to visually correlate high-risk radiomics zones with tumor regions like necrosis or cellular density.

## 5. Clinical Deployment & PACS Integration
*   **Objective**: Transition the Gradio demonstrator into a clinical decision support tool.
*   **Implementation**: Develop containerized APIs (using Docker) that can interface with clinical Picture Archiving and Communication Systems (PACS) via DICOM query/retrieve protocols, allowing automated GTV preprocessing and risk scoring.
