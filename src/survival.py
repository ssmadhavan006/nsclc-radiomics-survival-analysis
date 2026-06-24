import os
import time
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from lifelines import CoxPHFitter
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from typing import Dict, Any, List, Tuple

from src.utils import load_config
from src.model_utils import (
    plot_km_with_risk_table,
    bootstrap_cox_performance,
    calculate_time_dependent_roc,
    plot_calibration_curve_3yr,
    save_feature_meanings
)

logger = logging.getLogger("radiomics_pipeline")


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
        
        # Serialize model checkpoints (Task 4)
        scaler_selected = StandardScaler()
        scaler_selected.fit(df[selected_features])
        joblib.dump(scaler_selected, os.path.join(feat_dir, "scaler.joblib"))
        joblib.dump(cph_radiomics, os.path.join(feat_dir, "model_radiomics.joblib"))
        joblib.dump(cph_combined, os.path.join(feat_dir, "model_combined.joblib"))
        logger.info("Saved serialized model checkpoints to outputs/features/")

