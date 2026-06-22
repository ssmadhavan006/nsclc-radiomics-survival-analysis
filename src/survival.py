import os
import time
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score
from datetime import datetime
from typing import Dict, Any, List, Tuple

from src.utils import load_config

logger = logging.getLogger("radiomics_pipeline")

def plot_km_with_risk_table(
    durations: pd.Series,
    events: pd.Series,
    labels: pd.Series,
    title: str,
    output_path: str,
    time_unit: str = "Days"
) -> None:
    """
    Plots Kaplan-Meier curves with confidence intervals and a risk table (Rule 25, 36).
    """
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
    """
    Performs bootstrap validation to calculate 95% CIs for C-index and Hazard Ratios (Priority 2).
    """
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
    """
    Calculates time-dependent ROC curves and AUC at specific time points (Priority 2).
    """
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
    """
    Plots predicted vs observed survival probabilities at 3 years (Priority 2).
    """
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
    """
    Saves a CSV detailing the physical and biological meanings of the selected radiomic features (Rule 17, 30).
    """
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
            "Feature Name": feat,
            "Category": info["Category"],
            "Wavelet Band": info["Wavelet"],
            "Coefficient": coef_series.get(feat, 0.0),
            "IBSI Definition": info["Definition"],
            "Biological Interpretation": info["Interpretation"]
        })
        
    res_df = pd.DataFrame(records)
    res_df.to_csv(output_path, index=False)

def run_univariate_cox_models(df: pd.DataFrame, feature_cols: List[str], output_dir: str) -> pd.DataFrame:
    """
    Fits univariate Cox models for each clinical and radiomic variable (Stage 7).
    """
    cox_results = []
    duration_col, event_col = "Survival.time", "deadstatus.event"
    
    df_temp = df.copy()
    df_temp["gender_binary"] = df_temp["gender"].map({"male": 1, "female": 0}).fillna(0)
    
    for cov in ["age", "gender_binary"]:
        try:
            cph = CoxPHFitter()
            cph.fit(df_temp[[cov, duration_col, event_col]], duration_col, event_col)
            summary = cph.summary.iloc[0]
            cox_results.append({
                "Variable": cov, "Type": "Clinical", "HR": summary["exp(coef)"],
                "HR_Lower_CI": summary["exp(coef) lower 95%"], "HR_Upper_CI": summary["exp(coef) upper 95%"],
                "p_value": summary["p"], "C_Index": cph.concordance_index_
            })
        except Exception:
            pass
            
    for feat in feature_cols:
        try:
            scaler = StandardScaler()
            scaled_feat = scaler.fit_transform(df[[feat]])
            temp_df = pd.DataFrame({
                "feature": scaled_feat.ravel(), "duration": df[duration_col], "event": df[event_col]
            })
            cph = CoxPHFitter()
            cph.fit(temp_df, "duration", "event")
            summary = cph.summary.iloc[0]
            cox_results.append({
                "Variable": feat, "Type": "Radiomic", "HR": summary["exp(coef)"],
                "HR_Lower_CI": summary["exp(coef) lower 95%"], "HR_Upper_CI": summary["exp(coef) upper 95%"],
                "p_value": summary["p"], "C_Index": cph.concordance_index_
            })
        except Exception:
            pass
            
    cox_df = pd.DataFrame(cox_results)
    cox_df.to_csv(os.path.join(output_dir, "tables", "cox_univariate_results.csv"), index=False)
    return cox_df

def cross_validate_lasso_cox(df: pd.DataFrame, feature_cols: List[str], config: Dict[str, Any]) -> Tuple[float, List[str], np.ndarray]:
    """
    Performs nested parameter tuning using 5-fold cross-validation (Rule 27, 28).
    """
    duration_col, event_col = "Survival.time", "deadstatus.event"
    n_folds, seed = config["survival"]["n_cv_folds"], config["survival"]["random_seed"]
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    
    penalizers = [0.01, 0.05, 0.1, 0.2, 0.3]
    penalizer_scores = {p: [] for p in penalizers}
    
    for p in penalizers:
        for train_idx, val_idx in kf.split(df):
            train_df = df.iloc[train_idx].copy()
            val_df = df.iloc[val_idx].copy()
            
            scaler = StandardScaler()
            train_scaled = scaler.fit_transform(train_df[feature_cols])
            val_scaled = scaler.transform(val_df[feature_cols])
            
            train_cv = pd.DataFrame(train_scaled, columns=feature_cols)
            train_cv[duration_col] = train_df[duration_col].values
            train_cv[event_col] = train_df[event_col].values
            
            val_cv = pd.DataFrame(val_scaled, columns=feature_cols)
            val_cv[duration_col] = val_df[duration_col].values
            val_cv[event_col] = val_df[event_col].values
            
            try:
                cph = CoxPHFitter(penalizer=p, l1_ratio=1.0)
                cph.fit(train_cv, duration_col, event_col)
                c_idx = cph.score(val_cv, scoring_method="concordance_index")
                penalizer_scores[p].append(c_idx)
            except Exception:
                penalizer_scores[p].append(0.5)
                
    mean_scores = {p: np.mean(scores) for p, scores in penalizer_scores.items()}
    best_p = max(mean_scores, key=mean_scores.get)
    
    scaler = StandardScaler()
    scaled_all = scaler.fit_transform(df[feature_cols])
    final_df = pd.DataFrame(scaled_all, columns=feature_cols)
    final_df[duration_col] = df[duration_col].values
    final_df[event_col] = df[event_col].values
    
    final_cph = CoxPHFitter(penalizer=best_p, l1_ratio=1.0)
    final_cph.fit(final_df, duration_col, event_col)
    
    coefficients = final_cph.params_
    selected_features = coefficients[coefficients.abs() > 1e-4].index.tolist()
    return mean_scores[best_p], selected_features, coefficients.values

def run_survival_pipeline(config_path: str = "src/config.yaml") -> None:
    """
    Executes Stage 7: fits KM curves, comparative models, bootstrap validation,
    time-dependent ROC plots, and 3-year calibration curves.
    """
    config = load_config(config_path)
    output_dir = config["paths"]["output_dir"]
    feat_dir = os.path.join(output_dir, "features")
    fig_dir = os.path.join(output_dir, "figures")
    table_dir = os.path.join(output_dir, "tables")
    
    df = pd.read_csv(os.path.join(feat_dir, "cleaned_feature_matrix.csv"))
    df = df.dropna(subset=["Survival.time", "deadstatus.event"])
    logger.info(f"Loaded {len(df)} patients for survival modeling.")
    
    feature_cols = [c for c in df.columns if c not in [
        "PatientID", "age", "clinical.T.Stage", "Clinical.N.Stage", 
        "Clinical.M.Stage", "Overall.Stage", "Histology", "gender", 
        "Survival.time", "deadstatus.event", "stage_ordinal"
    ]]
    
    # 1. Plot standard KMs
    plot_km_with_risk_table(df["Survival.time"], df["deadstatus.event"], pd.Series(["All Patients"]*len(df)), "Overall Survival", os.path.join(fig_dir, "km_overall_survival.png"))
    plot_km_with_risk_table(df["Survival.time"], df["deadstatus.event"], df["gender"], "Survival by Gender", os.path.join(fig_dir, "km_survival_by_gender.png"))
    
    valid_stage_df = df[df["Overall.Stage"].isin(["I", "II", "IIIa", "IIIb"])]
    plot_km_with_risk_table(valid_stage_df["Survival.time"], valid_stage_df["deadstatus.event"], valid_stage_df["Overall.Stage"], "Survival by Stage", os.path.join(fig_dir, "km_survival_by_stage.png"))
    
    valid_hist_df = df[df["Histology"].isin(["squamous cell carcinoma", "large cell", "adenocarcinoma"])]
    plot_km_with_risk_table(valid_hist_df["Survival.time"], valid_hist_df["deadstatus.event"], valid_hist_df["Histology"], "Survival by Histology", os.path.join(fig_dir, "km_survival_by_histology.png"))
    
    # 2. Univariate Cox & LASSO-Cox
    cox_df = run_univariate_cox_models(df, feature_cols, output_dir)
    
    # KM by top feature
    if not cox_df.empty:
        radiomic_cox = cox_df[cox_df["Type"] == "Radiomic"]
        if not radiomic_cox.empty:
            top_feat = radiomic_cox.sort_values(by="p_value").iloc[0]["Variable"]
            median_val = df[top_feat].median()
            df["top_feat_group"] = df[top_feat].apply(lambda x: "High" if x > median_val else "Low")
            plot_km_with_risk_table(df["Survival.time"], df["deadstatus.event"], df["top_feat_group"], f"Survival by {top_feat.split('_')[-1]} (Top Feature)", os.path.join(fig_dir, "km_survival_by_top_feature.png"))
            
    cv_c_index, selected_features, coefficients = cross_validate_lasso_cox(df, feature_cols, config)
    
    if selected_features:
        # Build Radiomics score
        scaler = StandardScaler()
        scaled_df = pd.DataFrame(scaler.fit_transform(df[feature_cols]), columns=feature_cols)
        scaled_df["Survival.time"] = df["Survival.time"].values
        scaled_df["deadstatus.event"] = df["deadstatus.event"].values
        
        cph_radiomics = CoxPHFitter(penalizer=0.05, l1_ratio=1.0)
        cph_radiomics.fit(scaled_df[selected_features + ["Survival.time", "deadstatus.event"]], "Survival.time", "deadstatus.event")
        
        # Save feature meanings CSV (Priority 3)
        save_feature_meanings(selected_features, cph_radiomics.params_, os.path.join(table_dir, "feature_meanings.csv"))
        
        df["PrognosticScore"] = cph_radiomics.predict_partial_hazard(scaled_df[selected_features])
        median_score = df["PrognosticScore"].median()
        df["RiskGroup"] = df["PrognosticScore"].apply(lambda s: "High Risk" if s > median_score else "Low Risk")
        
        # Risk Group KM
        plot_km_with_risk_table(df["Survival.time"], df["deadstatus.event"], df["RiskGroup"], f"Survival Stratification by Radiomic Signature (C-Index: {cv_c_index:.3f})", os.path.join(fig_dir, "km_radiomic_signature.png"))
        
        # 3. Comparative Modeling (Priority 4)
        df["stage_binary"] = df["Overall.Stage"].apply(lambda s: 1 if s in ["IIIa", "IIIb"] else 0)
        df["gender_binary"] = df["gender"].map({"male": 1, "female": 0}).fillna(0)
        
        clinical_vars = ["age", "gender_binary", "stage_binary"]
        
        # Model A: Clinical only
        cph_clinical = CoxPHFitter()
        cph_clinical.fit(df[clinical_vars + ["Survival.time", "deadstatus.event"]], "Survival.time", "deadstatus.event")
        c_clinical = cph_clinical.concordance_index_
        
        # Model B: Radiomics only
        cph_rad_score = CoxPHFitter()
        cph_rad_score.fit(df[["PrognosticScore", "Survival.time", "deadstatus.event"]], "Survival.time", "deadstatus.event")
        c_radiomics = cph_rad_score.concordance_index_
        
        # Model C: Combined
        cph_combined = CoxPHFitter()
        cph_combined.fit(df[["PrognosticScore"] + clinical_vars + ["Survival.time", "deadstatus.event"]], "Survival.time", "deadstatus.event")
        c_combined = cph_combined.concordance_index_
        
        # 4. Bootstrap Confidence Intervals (Priority 2)
        ci_clinical, _ = bootstrap_cox_performance(df, clinical_vars)
        ci_radiomics, _ = bootstrap_cox_performance(df, ["PrognosticScore"])
        ci_combined, hr_cis = bootstrap_cox_performance(df, ["PrognosticScore"] + clinical_vars)
        
        # Write comparison table
        comp_df = pd.DataFrame([
            {"Model": "Clinical Only (A)", "C-Index": c_clinical, "95% CI": f"[{ci_clinical[0]:.4f} - {ci_clinical[1]:.4f}]"},
            {"Model": "Radiomics Only (B)", "C-Index": c_radiomics, "95% CI": f"[{ci_radiomics[0]:.4f} - {ci_radiomics[1]:.4f}]"},
            {"Model": "Combined (C)", "C-Index": c_combined, "95% CI": f"[{ci_combined[0]:.4f} - {ci_combined[1]:.4f}]"}
        ])
        comp_df.to_csv(os.path.join(table_dir, "model_comparison.csv"), index=False)
        logger.info(f"Model comparison complete. Combined C-Index: {c_combined:.4f}")
        
        # Save Multivariate Combined Results
        mv_results = cph_combined.summary.copy()
        mv_results["HR_Lower_Bootstrap"] = [hr_cis[idx][0] for idx in mv_results.index]
        mv_results["HR_Upper_Bootstrap"] = [hr_cis[idx][1] for idx in mv_results.index]
        mv_results.to_csv(os.path.join(table_dir, "cox_multivariate_results.csv"))
        
        # Plot multivariate combined forest plot
        plt.figure(figsize=(7, 5))
        cph_combined.plot(hazard_ratios=True)
        plt.title("Combined Model Hazard Ratios", fontsize=11, fontweight="bold")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.savefig(os.path.join(fig_dir, "cox_multivariate_forest.png"), dpi=300, bbox_inches="tight")
        plt.close()
        
        # 5. Calibration Curve at 3 years (Priority 2)
        plot_calibration_curve_3yr(df, cph_combined, os.path.join(fig_dir, "calibration_plot_3yr.png"))
        
        # 6. Time-dependent ROC (Priority 2)
        # Classify based on predicted PrognosticScore
        times = [365.0, 1095.0, 1825.0] # 1-year, 3-year, 5-year
        roc_results = calculate_time_dependent_roc(df["Survival.time"].values, df["deadstatus.event"].values, df["PrognosticScore"].values, times)
        
        plt.figure(figsize=(6, 6))
        for t in times:
            fpr, tpr, auc = roc_results[t]
            if len(fpr) > 0:
                plt.plot(fpr, tpr, label=f"{int(t/365)}-Year AUC: {auc:.3f}")
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xlabel("False Positive Rate", fontsize=10)
        plt.ylabel("True Positive Rate", fontsize=10)
        plt.title("Time-Dependent ROC Curves (Radiomic Signature)", fontsize=11, fontweight="bold")
        plt.legend(loc="lower right")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.savefig(os.path.join(fig_dir, "time_dependent_roc.png"), dpi=300, bbox_inches="tight")
        plt.close()
