# Radiomics Signature Feature Interpretation

This document details the 19 quantitative features selected by the regularized LASSO-Cox model to build the Prognostic Score signature.

| Programmatic Name | Friendly Name | Category | Coefficient | Wavelet | Physical & Biological Interpretation |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `original_shape_Elongation` | **Tumor Elongation** | Shape | $-0.0367$ | Original | Ratio of the minor axis to major axis of the tumor ellipsoid. Lower values indicate asymmetric, aggressive growth along anatomical structures. |
| `wavelet-LLH_firstorder_Maximum` | **Maximum Density (LLH)** | First-Order | $0.0711$ | LLH | Maximum voxel value after smoothing and horizontal high-pass filtering. Reflects localized high-density active components. |
| `wavelet-LLH_firstorder_Median` | **Median Density (LLH)** | First-Order | $0.0543$ | LLH | Median voxel value after LLH filtering. Represents smoothed central tumor density. |
| `wavelet-LLH_firstorder_TotalEnergy` | **Total Density Energy (LLH)** | First-Order | $0.0489$ | LLH | Volume-scaled sum of squared voxel values. Combines physical tumor volume and overall density burden. |
| `wavelet-LHL_firstorder_Skewness` | **Density Skewness (LHL)** | First-Order | $-0.0152$ | LHL | Asymmetry of intensity distributions. Indicates localized areas of necrosis or viability. |
| `wavelet-LHL_glrlm_LongRunHighGrayLevelEmphasis` | **Long Dense Runs Emphasis (LHL)** | GLRLM | $0.0321$ | LHL | Emphasizes long runs of high-intensity voxels. Represents large, continuous zones of active solid tumor mass. |
| `wavelet-LHL_gldm_DependenceVariance` | **Texture Dependency Variance (LHL)** | GLDM | $0.0211$ | LHL | Variance of local voxel dependencies. Reflects structural complexity and microenvironmental heterogeneity. |
| `wavelet-LHL_ngtdm_Strength` | **Texture Strength (LHL)** | NGTDM | $0.0189$ | LHL | Measures contrast transitions. High values represent distinct interfaces between different tumor tissue types. |
| `wavelet-LHH_firstorder_Kurtosis` | **Density Peakedness (LHH)** | First-Order | $0.0125$ | LHH | Peakedness of intensity distributions in the high-frequency band. Highlights small spots of necrosis or calcification. |
| `wavelet-LHH_firstorder_Maximum` | **Maximum Density (LHH)** | First-Order | $0.0432$ | LHH | Localized maximum density peak under horizontal/vertical high-frequency bands. |
| `wavelet-HHL_firstorder_Skewness` | **Density Skewness (HHL)** | First-Order | $-0.0112$ | HHL | Asymmetry of high-frequency intensity distributions. Reflects directional asymmetry in cellular density. |
| `wavelet-HHL_glszm_SizeZoneNonUniformity` | **Size Zone Non-Uniformity (HHL)** | GLSZM | $0.0892$ | HHL | Measures size zone volume variability. High values indicate highly fragmented and heterogeneous tumor zones. |
| `wavelet-HHH_glszm_GrayLevelNonUniformity` | **Zone Intensity Heterogeneity (HHH)** | GLSZM | $0.0612$ | HHH | Intensity variability across size zones. Indicates spatial heterogeneity in active cellular density. |
| `wavelet-HHH_glszm_SizeZoneNonUniformity` | **Size Zone Non-Uniformity (HHH)** | GLSZM | $0.0754$ | HHH | Size zone volume variability under high-frequency HHH filtering. Represents structural fragmentation. |
| `wavelet-HHH_glszm_SmallAreaHighGrayLevelEmphasis` | **Small Dense Zones Emphasis (HHH)** | GLSZM | $0.1616$ | HHH | Emphasizes small, high-density voxel clusters. Reflects active microcalcifications or highly proliferative tumor sub-regions. |
| `wavelet-HHH_glszm_ZoneVariance` | **Zone Size Variance (HHH)** | GLSZM | $0.0543$ | HHH | Variance of zone volumes. Indicates complex structural patterns mixing large and small density zones. |
| `wavelet-LLL_firstorder_Minimum` | **Minimum Density (LLL)** | First-Order | $0.0212$ | LLL | Minimum voxel intensity after 3D low-pass smoothing. Typically corresponds to necrosis or fluid inside the tumor core. |
| `wavelet-LLL_firstorder_Range` | **Density Range (LLL)** | First-Order | $-0.0948$ | LLL | Intensity range in low-pass smoothed volumes. Represents broad density gradients across necrotic and solid tumor components. |
| `wavelet-LLL_gldm_SmallDependenceHighGrayLevelEmphasis` | **Small Dense Dependencies (LLL)** | GLDM | $0.0452$ | LLL | Emphasizes small, high-density dependencies. Indicates fine-grained cellular clusters inside smoothed volumes. |
