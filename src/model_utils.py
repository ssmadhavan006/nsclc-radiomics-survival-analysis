import os
import yaml
import logging
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
from datetime import datetime
from typing import Dict, Any, List, Tuple
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score

logger = logging.getLogger("radiomics_pipeline")

FEATURE_FRIENDLY_NAMES = {
    "original_shape_Elongation": "Tumor Elongation",
    "wavelet-LLH_firstorder_Maximum": "Maximum Density (LLH)",
    "wavelet-LLH_firstorder_Median": "Median Density (LLH)",
    "wavelet-LLH_firstorder_TotalEnergy": "Total Density Energy (LLH)",
    "wavelet-LHL_firstorder_Skewness": "Density Skewness (LHL)",
    "wavelet-LHL_glrlm_LongRunHighGrayLevelEmphasis": "Long Dense Runs Emphasis (LHL)",
    "wavelet-LHL_gldm_DependenceVariance": "Texture Dependency Variance (LHL)",
    "wavelet-LHL_ngtdm_Strength": "Texture Strength (LHL)",
    "wavelet-LHH_firstorder_Kurtosis": "Density Peakedness (LHH)",
    "wavelet-LHH_firstorder_Maximum": "Maximum Density (LHH)",
    "wavelet-HHL_firstorder_Skewness": "Density Skewness (HHL)",
    "wavelet-HHL_glszm_SizeZoneNonUniformity": "Size Zone Non-Uniformity (HHL)",
    "wavelet-HHH_glszm_GrayLevelNonUniformity": "Zone Intensity Heterogeneity (HHH)",
    "wavelet-HHH_glszm_SizeZoneNonUniformity": "Size Zone Non-Uniformity (HHH)",
    "wavelet-HHH_glszm_SmallAreaHighGrayLevelEmphasis": "Small Dense Zones Emphasis (HHH)",
    "wavelet-HHH_glszm_ZoneVariance": "Zone Size Variance (HHH)",
    "wavelet-LLL_firstorder_Minimum": "Minimum Density (LLL)",
    "wavelet-LLL_firstorder_Range": "Density Range (LLL)",
    "wavelet-LLL_gldm_SmallDependenceHighGrayLevelEmphasis": "Small Dense Dependencies (LLL)"
}

def get_friendly_feature_name(feat_name: str) -> str:
    """
    Returns a clean, human-readable name for a programmatic radiomics feature name.
    """
    if feat_name in FEATURE_FRIENDLY_NAMES:
        return FEATURE_FRIENDLY_NAMES[feat_name]
    
    parts = feat_name.split("_")
    if len(parts) >= 3:
        metric = parts[-1]
        category = parts[-2]
        wavelet = parts[-3].replace("wavelet-", "")
        return f"{metric} ({category} {wavelet})"
    elif len(parts) == 2:
        return f"{parts[1]} ({parts[0]})"
    return feat_name

def format_ordinal(n: float) -> str:
    """Formats a percentage rank to its ordinal string representation."""
    val = int(round(n))
    if 11 <= (val % 100) <= 13:
        suffix = 'th'
    else:
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(val % 10, 'th')
    return f"{val}{suffix} percentile"

def format_risk_markdown(risk_cat: str, hazard_str: str, percentile_str: str) -> str:
    """Generates styled HTML representation of risk stratifications and confidence metrics."""
    color = "#f87171" if risk_cat == "High Risk" else "#4ade80"
    return f"""
<div style='background-color: #0f172a; padding: 14px; border-radius: 8px; border: 1px solid #334155;'>
    <div style='margin-bottom: 12px;'>
        <span style='color: #94a3b8; font-size: 11px; text-transform: uppercase; font-weight: bold;'>Risk Category</span><br/>
        <strong style='color: {color}; font-size: 20px;'>{risk_cat}</strong>
    </div>
    <div style='margin-bottom: 12px; border-top: 1px solid #334155; padding-top: 8px;'>
        <span style='color: #94a3b8; font-size: 11px; text-transform: uppercase; font-weight: bold;'>Predicted Hazard Deviation</span><br/>
        <strong style='color: #f1f5f9; font-size: 14px;'>{hazard_str}</strong>
    </div>
    <div style='border-top: 1px solid #334155; padding-top: 8px;'>
        <span style='color: #94a3b8; font-size: 11px; text-transform: uppercase; font-weight: bold;'>Cohort Percentile Rank</span><br/>
        <strong style='color: #38bdf8; font-size: 14px;'>{percentile_str}</strong>
    </div>
</div>
"""

def plot_km_with_risk_table(
    durations: pd.Series,
    events: pd.Series,
    labels: pd.Series,
    title: str,
    output_path: str,
    time_unit: str = "Days"
) -> None:
    """Plots Kaplan-Meier curves with confidence intervals and a risk table."""
    unique_groups = sorted(labels.unique())
    fig = plt.figure(figsize=(10, 8))
    gs = plt.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.3)
    ax_km = fig.add_subplot(gs[0])
    ax_table = fig.add_subplot(gs[1])
    
    kmfs = {}
    colors = sns.color_palette("Set1", len(unique_groups))
    
    for idx, group in enumerate(unique_groups):
        mask = (labels == group)
        kmf = KaplanMeierFitter()
        kmf.fit(durations[mask], events[mask], label=str(group))
        kmf.plot_survival_function(
            ax=ax_km, color=colors[idx], show_censors=True, ci_show=True,
            censor_styles={"ms": 6, "mew": 1.5}
        )
        kmfs[group] = kmf
        
    ax_km.set_title(title, fontsize=12, fontweight="bold")
    ax_km.set_xlabel("")
    ax_km.set_ylabel("Survival Probability", fontsize=10)
    ax_km.grid(True, linestyle="--", alpha=0.5)
    
    if len(unique_groups) == 2:
        m1, m2 = (labels == unique_groups[0]), (labels == unique_groups[1])
        res = logrank_test(durations[m1], durations[m2], event_observed_A=events[m1], event_observed_B=events[m2])
        ax_km.text(
            0.05, 0.05, f"Log-rank p: {res.p_value:.3e}" if res.p_value < 0.001 else f"Log-rank p: {res.p_value:.4f}",
            transform=ax_km.transAxes, fontsize=9, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
        )
        
    time_ticks = np.linspace(0, durations.max(), 6)
    ax_km.set_xticks(time_ticks)
    ax_km.set_xticklabels([f"{int(t)}" for t in time_ticks])
    
    ax_table.axis('off')
    table_data = []
    row_labels = [f"Group {g}" for g in unique_groups]
    
    for group in unique_groups:
        row_counts = []
        for tick in time_ticks:
            at_risk = kmfs[group].event_table.loc[tick:]["at_risk"].iloc[0] if tick in kmfs[group].event_table.index else 0
            if tick > 0 and at_risk == 0:
                past_events = kmfs[group].event_table.loc[:tick]
                at_risk = past_events["at_risk"].iloc[-1] - past_events["removed"].iloc[-1] if not past_events.empty else 0
            row_counts.append(int(at_risk))
        table_data.append(row_counts)
        
    table = ax_table.table(
        cellText=table_data, rowLabels=row_labels, colLabels=[f"{int(t)} {time_unit[:3]}" for t in time_ticks],
        cellLoc='center', loc='center'
    )
    table.scale(1, 1.3)
    table.set_fontsize(8)
    
    fig.text(0.02, 0.02, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}", fontsize=7, style="italic")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

def bootstrap_cox_performance(
    df: pd.DataFrame,
    covariates: List[str],
    n_bootstrap: int = 1000,
    seed: int = 42
) -> Tuple[Tuple[float, float], Dict[str, Tuple[float, float]]]:
    """Performs bootstrap validation to calculate 95% CIs for C-index and Hazard Ratios."""
    np.random.seed(seed)
    c_indexes = []
    hr_samples = {cov: [] for cov in covariates}
    n_samples = len(df)
    
    temp_df = df[covariates + ["Survival.time", "deadstatus.event"]].copy()
    
    for _ in range(n_bootstrap):
        boot_df = temp_df.sample(n=n_samples, replace=True).reset_index(drop=True)
        if boot_df["deadstatus.event"].sum() in [0, len(boot_df)]:
            continue
            
        cph = CoxPHFitter()
        try:
            cph.fit(boot_df, "Survival.time", "deadstatus.event")
            c_indexes.append(cph.concordance_index_)
            for cov in covariates:
                hr_samples[cov].append(np.exp(cph.params_[cov]))
        except Exception:
            continue
            
    c_index_ci = (np.percentile(c_indexes, 2.5), np.percentile(c_indexes, 97.5))
    hr_cis = {cov: (np.percentile(hr_samples[cov], 2.5), np.percentile(hr_samples[cov], 97.5)) for cov in covariates if hr_samples[cov]}
    return c_index_ci, hr_cis

def calculate_time_dependent_roc(
    durations: np.ndarray,
    events: np.ndarray,
    predictions: np.ndarray,
    time_points: List[float]
) -> Dict[float, Tuple[np.ndarray, np.ndarray, float]]:
    """Calculates time-dependent ROC curves and AUC at specific time points."""
    results = {}
    for t in time_points:
        labels, scores = [], []
        for d, e, s in zip(durations, events, predictions):
            if d > t:
                labels.append(0)
                scores.append(s)
            elif d <= t and e == 1:
                labels.append(1)
                scores.append(s)
        if len(np.unique(labels)) > 1:
            fpr, tpr, _ = roc_curve(labels, scores)
            auc = roc_auc_score(labels, scores)
            results[t] = (fpr, tpr, auc)
        else:
            results[t] = (np.array([]), np.array([]), np.nan)
    return results

def plot_calibration_curve_3yr(
    df: pd.DataFrame,
    cph: CoxPHFitter,
    output_path: str,
    time_point: float = 1095,
    n_bins: int = 5
) -> None:
    """Plots predicted vs observed survival probabilities at 3 years."""
    surv_probs = cph.predict_survival_function(df, times=[time_point]).T.iloc[:, 0]
    df_cal = df.copy()
    df_cal["pred_surv"] = surv_probs
    df_cal["bin"] = pd.qcut(df_cal["pred_surv"], q=n_bins, labels=False, duplicates="drop")
    
    pred_vals, obs_vals, obs_cis = [], [], []
    for b in sorted(df_cal["bin"].unique()):
        bin_df = df_cal[df_cal["bin"] == b]
        pred_vals.append(bin_df["pred_surv"].mean())
        
        kmf = KaplanMeierFitter()
        kmf.fit(bin_df["Survival.time"], bin_df["deadstatus.event"])
        obs_mean = kmf.survival_function_at_times(time_point).values[0]
        obs_vals.append(obs_mean)
        
        ci = kmf.confidence_interval_
        closest_idx = np.abs(ci.index.values - time_point).argmin()
        obs_cis.append((ci.iloc[closest_idx].iloc[0], ci.iloc[closest_idx].iloc[1]))
        
    plt.figure(figsize=(6, 6))
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Ideal (Perfect)")
    yerr = np.array([[obs - ci[0], ci[1] - obs] for obs, ci in zip(obs_vals, obs_cis)]).T
    plt.errorbar(pred_vals, obs_vals, yerr=yerr, fmt="o-", color="navy", label="Model Calibration", elinewidth=1.5, capsize=3)
    
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel("Predicted 3-Year Survival Probability", fontsize=10)
    plt.ylabel("Observed 3-Year Survival Probability", fontsize=10)
    plt.title("Model Calibration Plot at 3 Years", fontsize=12, fontweight="bold")
    plt.legend(loc="upper left")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

def save_feature_meanings(selected_features: List[str], coef_series: pd.Series, output_path: str) -> None:
    """Saves a CSV detailing the physical and biological meanings of the selected features."""
    meanings = {
        "original_shape_Elongation": {
            "Category": "Shape",
            "Wavelet": "Original",
            "Definition": "Ratio of the minor axis to the major axis of the tumor ellipsoid.",
            "Interpretation": "Measures elongation. Lower values indicate highly elongated tumors, reflecting asymmetric, aggressive growth along anatomical structures."
        },
        "wavelet-LLH_firstorder_Maximum": {
            "Category": "First-Order",
            "Wavelet": "LLH",
            "Definition": "Maximum voxel value after LLH wavelet filtering.",
            "Interpretation": "Reflects localized high-density peaks (e.g. solid components/calcification) in a low-frequency background."
        },
        "wavelet-LLH_firstorder_Median": {
            "Category": "First-Order",
            "Wavelet": "LLH",
            "Definition": "Median voxel value after LLH wavelet filtering.",
            "Interpretation": "Represents the central density of the tumor after smoothing high-frequency noise horizontally/vertically."
        },
        "wavelet-LLH_firstorder_TotalEnergy": {
            "Category": "First-Order",
            "Wavelet": "LLH",
            "Definition": "Sum of squared voxel values scaled by voxel volume.",
            "Interpretation": "Combines tumor volume and density. Large, dense tumors produce high values, indicating overall tumor burden."
        },
        "wavelet-LHL_firstorder_Skewness": {
            "Category": "First-Order",
            "Wavelet": "LHL",
            "Definition": "Asymmetry of voxel intensity distribution about its mean.",
            "Interpretation": "Indicates spatial asymmetry in tissue composition (e.g. localized necrotic vs viable tumor regions)."
        },
        "wavelet-LHL_glrlm_LongRunHighGrayLevelEmphasis": {
            "Category": "GLRLM",
            "Wavelet": "LHL",
            "Definition": "Joint emphasis of long runs and high gray-level values.",
            "Interpretation": "Represents large, continuous zones of high density (e.g. active solid tumor mass)."
        },
        "wavelet-LHL_gldm_DependenceVariance": {
            "Category": "GLDM",
            "Wavelet": "LHL",
            "Definition": "Variance of local voxel dependency counts.",
            "Interpretation": "Indicates variation in local texture density, suggesting micro-environmental structural heterogeneity."
        },
        "wavelet-LHL_ngtdm_Strength": {
            "Category": "NGTDM",
            "Wavelet": "LHL",
            "Definition": "Measures texture contrast and sharpness.",
            "Interpretation": "Reflects clear interfaces and high contrast transitions between different tumor tissue types."
        },
        "wavelet-LHH_firstorder_Kurtosis": {
            "Category": "First-Order",
            "Wavelet": "LHH",
            "Definition": "Peakedness of the intensity distribution.",
            "Interpretation": "High values indicate intensity variations are dominated by extreme values (e.g. microcalcifications/necrotic spots)."
        },
        "wavelet-LHH_firstorder_Maximum": {
            "Category": "First-Order",
            "Wavelet": "LHH",
            "Definition": "Maximum voxel value in the horizontal high-frequency band.",
            "Interpretation": "Reflects fine-grained high-intensity details or sharp horizontal density transitions."
        },
        "wavelet-HHL_firstorder_Skewness": {
            "Category": "First-Order",
            "Wavelet": "HHL",
            "Definition": "Intensity asymmetry in the HHL wavelet band.",
            "Interpretation": "Reflects directional asymmetry in cellular density or vascular patterns."
        },
        "wavelet-HHL_glszm_SizeZoneNonUniformity": {
            "Category": "GLSZM",
            "Wavelet": "HHL",
            "Definition": "Variability of size zone volumes throughout the tumor.",
            "Interpretation": "High values indicate a wide variety of zone sizes, representing highly fragmented, heterogeneous tumor zones."
        },
        "wavelet-HHH_glszm_GrayLevelNonUniformity": {
            "Category": "GLSZM",
            "Wavelet": "HHH",
            "Definition": "Variability of gray-level intensities across size zones.",
            "Interpretation": "Indicates significant variation in density levels across zones, reflecting intratumoral density heterogeneity."
        },
        "wavelet-HHH_glszm_SizeZoneNonUniformity": {
            "Category": "GLSZM",
            "Wavelet": "HHH",
            "Definition": "Variability of size zone volumes under HHH filtering.",
            "Interpretation": "Highly heterogeneous distribution of fine-texture zone sizes, indicating structural fragmentation."
        },
        "wavelet-HHH_glszm_SmallAreaHighGrayLevelEmphasis": {
            "Category": "GLSZM",
            "Wavelet": "HHH",
            "Definition": "Emphasizes small zones of high gray-level values.",
            "Interpretation": "Reflects small, highly dense clusters of active cells or microcalcifications, indicating active proliferation."
        },
        "wavelet-HHH_glszm_ZoneVariance": {
            "Category": "GLSZM",
            "Wavelet": "HHH",
            "Definition": "Variance in zone size volumes.",
            "Interpretation": "Indicates complex structural patterns with a mix of very small and very large density zones."
        },
        "wavelet-LLL_firstorder_Minimum": {
            "Category": "First-Order",
            "Wavelet": "LLL",
            "Definition": "Minimum voxel value in the smoothed LLL band.",
            "Interpretation": "Represents the lowest density in the smoothed volume, typically corresponding to necrosis or fluid."
        },
        "wavelet-LLL_firstorder_Range": {
            "Category": "First-Order",
            "Wavelet": "LLL",
            "Definition": "Difference between max and min voxel values in the LLL band.",
            "Interpretation": "Measures the overall range of macro-intensities, reflecting broad density gradients."
        },
        "wavelet-LLL_gldm_SmallDependenceHighGrayLevelEmphasis": {
            "Category": "GLDM",
            "Wavelet": "LLL",
            "Definition": "Emphasizes small dependencies with high gray-level values.",
            "Interpretation": "Indicates fine-grained high-density textures in the smoothed volume, reflecting clusters of viable tumor cells."
        }
    }
    
    records = []
    for feat in selected_features:
        info = meanings.get(feat, {
            "Category": "Unknown",
            "Wavelet": "Unknown",
            "Definition": "Custom texture/intensity metric.",
            "Interpretation": "Reflects spatial density variation."
        })
        records.append({
            "Feature Name": get_friendly_feature_name(feat),
            "Category": info["Category"],
            "Wavelet Band": info["Wavelet"],
            "Coefficient": coef_series.get(feat, 0.0),
            "IBSI Definition": info["Definition"],
            "Biological Interpretation": info["Interpretation"]
        })
        
    res_df = pd.DataFrame(records)
    res_df.to_csv(output_path, index=False)

def generate_roi_image(patient_id: str, ct_dir: str = None, seg_path: str = None) -> str:
    """Generates side-by-side axial crop of CT and CT+GTV Overlay with contour outline."""
    cache_dir = "outputs/figures/cache"
    os.makedirs(cache_dir, exist_ok=True)
    out_path = os.path.join(cache_dir, f"roi_{patient_id}.png")
    
    # Check if cached image already exists
    if os.path.exists(out_path) and ct_dir is None:
        return out_path
        
    try:
        # Load config
        with open("src/config.yaml", "r") as f:
            config = yaml.safe_load(f)
            
        # Discover paths for demo patient if not provided
        if ct_dir is None or seg_path is None:
            patient_dir = f"dataset/NSCLC-Radiomics/{patient_id}"
            from src.data_ingestion import find_patient_series
            ct_dir, seg_path, _ = find_patient_series(patient_dir)
            if ct_dir is None or seg_path is None:
                logger.error(f"Failed to find DICOM files for patient {patient_id}")
                return None
                
        # Run preprocessing to get coordinate-aligned isotropic CT and Mask
        from src.preprocessing import preprocess_case
        preprocessed_ct, preprocessed_mask = preprocess_case(ct_dir, seg_path, config)
        
        # Get numpy arrays
        ct_arr = sitk.GetArrayFromImage(preprocessed_ct)
        mask_arr = sitk.GetArrayFromImage(preprocessed_mask)
        
        # Find slice with maximum tumor area
        slice_sums = np.sum(mask_arr, axis=(1, 2))
        max_slice_idx = int(np.argmax(slice_sums))
        
        if slice_sums[max_slice_idx] == 0:
            max_slice_idx = ct_arr.shape[0] // 2
            
        ct_slice = ct_arr[max_slice_idx]
        mask_slice = mask_arr[max_slice_idx]
        
        # Normalize CT slice for display (-1000 to 400 HU)
        min_hu, max_hu = -1000.0, 400.0
        ct_norm = np.clip(ct_slice, min_hu, max_hu)
        ct_norm = ((ct_norm - min_hu) / (max_hu - min_hu) * 255.0).astype(np.uint8)
        
        # Create RGB image
        rgb_slice = np.stack([ct_norm, ct_norm, ct_norm], axis=-1)
        
        # Crop region centered around the GTV centroid on this slice
        rows, cols = np.where(mask_slice > 0)
        if len(rows) > 0:
            center_y = int(np.mean(rows))
            center_x = int(np.mean(cols))
            
            crop_size = 160
            y_min = max(0, center_y - crop_size // 2)
            y_max = min(rgb_slice.shape[0], center_y + crop_size // 2)
            x_min = max(0, center_x - crop_size // 2)
            x_max = min(rgb_slice.shape[1], center_x + crop_size // 2)
            
            # Adjust if crop size is smaller due to boundaries
            if y_max - y_min < crop_size:
                if y_min == 0:
                    y_max = min(rgb_slice.shape[0], crop_size)
                else:
                    y_min = max(0, rgb_slice.shape[0] - crop_size)
            if x_max - x_min < crop_size:
                if x_min == 0:
                    x_max = min(rgb_slice.shape[1], crop_size)
                else:
                    x_min = max(0, rgb_slice.shape[1] - crop_size)
        else:
            center_y, center_x = rgb_slice.shape[0] // 2, rgb_slice.shape[1] // 2
            crop_size = 160
            y_min, y_max = center_y - crop_size // 2, center_y + crop_size // 2
            x_min, x_max = center_x - crop_size // 2, center_x + crop_size // 2
            
        # Draw GTV overlay on a copy of the RGB image
        rgb_overlay = rgb_slice.copy()
        mask_mask = (mask_slice > 0)
        alpha = 0.35
        overlay_color = [255, 0, 0] # Red
        rgb_overlay[mask_mask] = (alpha * np.array(overlay_color) + (1 - alpha) * rgb_slice[mask_mask]).astype(np.uint8)
        
        # Crop both original and overlay
        cropped_orig = rgb_slice[y_min:y_max, x_min:x_max]
        cropped_overlay = rgb_overlay[y_min:y_max, x_min:x_max]
        
        # Add a sharp boundary outline to the tumor in the overlay image
        edge_mask = mask_slice ^ ndimage.binary_erosion(mask_slice)
        cropped_edge = edge_mask[y_min:y_max, x_min:x_max]
        cropped_overlay[cropped_edge > 0] = [255, 50, 50] # Bright red outline
        
        # Plot and save
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        ax1.imshow(cropped_orig)
        ax1.set_title("Original CT (Axial Zoom)", fontsize=12, fontweight="bold", color="#f1f5f9")
        ax1.axis("off")
        
        ax2.imshow(cropped_overlay)
        ax2.set_title("Tumor Mask Overlay (GTV)", fontsize=12, fontweight="bold", color="#f1f5f9")
        ax2.axis("off")
        
        # Set figure background to match Gradio dark theme
        fig.patch.set_facecolor("#1e293b")
        plt.tight_layout()
        
        fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor(), edgecolor='none')
        plt.close(fig)
        
        return out_path
    except Exception as e:
        logger.error(f"Failed to generate GTV overlay: {str(e)}")
        return None

def process_patient_directory(ct_dir, seg_path, selected_features):
    """Processes a patient CT folder and SEG file to extract the features on-the-fly."""
    try:
        # Load master config
        with open("src/config.yaml", "r") as f:
            config = yaml.safe_load(f)
            
        # Run ingestion & preprocessing steps
        from src.preprocessing import preprocess_case
        preprocessed_ct, preprocessed_mask = preprocess_case(ct_dir, seg_path, config)
        
        # Check GTV voxel overlaps
        mask_arr = sitk.GetArrayFromImage(preprocessed_mask)
        voxel_count = np.sum(mask_arr)
        if voxel_count < 50:
            return f"Error: Resampled GTV mask only has {voxel_count} voxels. Must contain at least 50.", None
            
        # Run PyRadiomics feature extractor
        from src.feature_extraction import get_radiomics_extractor
        extractor = get_radiomics_extractor(config)
        feature_vector = extractor.execute(preprocessed_ct, preprocessed_mask)
        
        # Filter out features
        extracted_features = {}
        for feat in selected_features:
            val = feature_vector.get(feat, 0.0)
            if hasattr(val, "item"):
                extracted_features[feat] = float(val.item())
            else:
                extracted_features[feat] = float(val)
                
        return "Success", extracted_features
        
    except Exception as e:
        return f"Pipeline execution failed: {str(e)}", None

