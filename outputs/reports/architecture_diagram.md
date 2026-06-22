# Radiomics Pipeline Architecture Diagram

This document contains a Mermaid flowchart depicting the high-level pipeline architecture from clinical and image data ingestion to prognostic validation.

```mermaid
graph TD
    subgraph Ingestion ["1. Data Ingestion & QC"]
        A[Clinical Metadata CSV] --> C[Ingestion Manifest Builder]
        B[CT/SEG DICOM Folders] --> C
        C --> D{Coordinate QC}
        D -->|Fail| E[Exclude: failed_cases.csv]
        D -->|Pass| F[Valid Cohort Manifest]
    end
    subgraph Preprocessing ["2. Image Preprocessing"]
        F --> G[HU Intensity Clipping: -1000 to 400 HU]
        G --> H[Isotropic Resampling: 1.0mm^3 B-Spline / NN]
        H --> I{Spatial Overlap Verification}
        I -->|Fail| E
        I -->|Pass| J[Coordinate-Aligned CT & Mask]
    end
    subgraph Extraction ["3. Feature Extraction"]
        J --> K[Parallel PyRadiomics Loop: LokyBackend]
        K --> L[Raw Features: raw_features_all_patients.csv]
    end
    subgraph Engineering ["4. Feature Engineering"]
        L --> M[Variance Filter: threshold = 0.01]
        M --> N[Spearman Redundancy Filter: threshold = 0.95]
        N --> O[Clean Features: cleaned_feature_matrix.csv]
    end
    subgraph Analysis ["5. Statistical Analysis"]
        O --> P[Univariate Association Tests]
        P --> Q[Benjamini-Hochberg FDR Correction]
        Q --> R[Stage-Associated Phenotypes]
    end
    subgraph Prognostic ["6. Survival Modeling"]
        O --> S[5-Fold Cross-Validated LASSO-Cox]
        S --> T[19 Selected Signature Features]
        T --> U[PrognosticScore Signature]
    end
    subgraph Validation ["7. Validation & Evaluation"]
        U --> V[Kaplan-Meier Risk Stratification]
        U --> W[1000x Bootstrap CIs & Multivariate Cox]
        U --> X[Time-Dependent ROC AUCs: 1, 3, 5-Year]
        U --> Y[3-Year Calibration Quintiles]
    end

    style Ingestion fill:#f9f,stroke:#333,stroke-width:2px
    style Preprocessing fill:#bbf,stroke:#333,stroke-width:2px
    style Extraction fill:#bfb,stroke:#333,stroke-width:2px
    style Engineering fill:#ffb,stroke:#333,stroke-width:2px
    style Analysis fill:#fbb,stroke:#333,stroke-width:2px
    style Prognostic fill:#fbf,stroke:#333,stroke-width:2px
    style Validation fill:#bff,stroke:#333,stroke-width:2px
```
