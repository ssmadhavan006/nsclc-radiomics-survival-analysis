import os
import time
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from datetime import datetime
from typing import Dict, Any, Tuple, List

from src.utils import load_config

logger = logging.getLogger("radiomics_pipeline")

def clean_and_prepare_features(
    raw_features_path: str,
    clinical_csv_path: str,
    config: Dict[str, Any]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Cleans raw radiomics features, merges with clinical data, handles missing values,
    and applies variance and correlation filtering (Stage 5).
    
    Args:
        raw_features_path: Path to the raw features CSV.
        clinical_csv_path: Path to the clinical CSV metadata.
        config: Pipeline configuration dictionary.
        
    Returns:
        Tuple: (cleaned_features_merged_df, dropped_features_df)
    """
    logger.info("Starting Stage 5: Data Cleaning & Feature Engineering")
    
    # 1. Load data
    raw_df = pd.read_csv(raw_features_path)
    clinical_df = pd.read_csv(clinical_csv_path)
    
    # Remove diagnostics columns
    diag_cols = [col for col in raw_df.columns if col.startswith("diagnostics_")]
    features_df = raw_df.drop(columns=diag_cols)
    
    # Merge datasets
    merged_df = pd.merge(features_df, clinical_df, on="PatientID", how="inner")
    logger.info(f"Merged raw features and clinical data. Matched patients: {len(merged_df)}")
    
    # 2. Handle missing values
    # Feature imputation (median)
    feature_cols = [col for col in features_df.columns if col != "PatientID"]
    for col in feature_cols:
        if merged_df[col].isnull().any():
            median_val = merged_df[col].median()
            merged_df[col] = merged_df[col].fillna(median_val)
            
    # Clinical variable imputation/handling
    merged_df["Histology"] = merged_df["Histology"].fillna("Unknown").astype(str).str.strip()
    merged_df["Overall.Stage"] = merged_df["Overall.Stage"].fillna("Unknown").astype(str).str.strip()
    
    if merged_df["age"].isnull().any():
        median_age = merged_df["age"].median()
        merged_df["age"] = merged_df["age"].fillna(median_age)
        
    # Save clinical summary table
    table_dir = os.path.join(config["paths"]["output_dir"], "tables")
    os.makedirs(table_dir, exist_ok=True)
    
    clinical_cols = ["age", "clinical.T.Stage", "Clinical.N.Stage", "Clinical.M.Stage", "Overall.Stage", "Histology", "gender"]
    clinical_desc = merged_df[clinical_cols].describe(include='all')
    clinical_desc.to_csv(os.path.join(table_dir, "clinical_summary.csv"))
    
    # 3. Variance filtering (Rule 14)
    var_threshold = config["analysis"]["variance_threshold"]
    dropped_log = []
    
    variances = merged_df[feature_cols].var()
    low_var_cols = variances[variances < var_threshold].index.tolist()
    
    for col in low_var_cols:
        dropped_log.append({
            "Feature": col,
            "Reason": "Low Variance",
            "Value": variances[col]
        })
        
    remaining_features = [col for col in feature_cols if col not in low_var_cols]
    logger.info(f"Removed {len(low_var_cols)} low-variance features (threshold {var_threshold}).")
    
    # 4. Correlation filtering (Rule 14)
    corr_threshold = config["analysis"]["correlation_threshold"]
    corr_matrix = merged_df[remaining_features].corr(method="spearman").abs()
    
    upper_tri = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    to_drop = set()
    
    # Average absolute correlation of each feature to all others
    mean_abs_corr = corr_matrix.mean()
    
    for col in upper_tri.columns:
        # Find features highly correlated with 'col'
        high_corr_features = upper_tri.index[upper_tri[col] > corr_threshold].tolist()
        if high_corr_features:
            # Pairwise redundancy: compare average correlation of col and the others
            # Keep the one with lower average correlation to others
            for f in high_corr_features:
                if mean_abs_corr[col] > mean_abs_corr[f]:
                    to_drop.add(col)
                else:
                    to_drop.add(f)
                    
    for col in to_drop:
        dropped_log.append({
            "Feature": col,
            "Reason": "High Correlation",
            "Value": mean_abs_corr[col]
        })
        
    final_features = [col for col in remaining_features if col not in to_drop]
    logger.info(f"Removed {len(to_drop)} highly-correlated features (threshold {corr_threshold}).")
    logger.info(f"Final feature count: {len(final_features)}")
    
    # Save dropped features log
    feat_dir = os.path.join(config["paths"]["output_dir"], "features")
    dropped_df = pd.DataFrame(dropped_log)
    dropped_df.to_csv(os.path.join(feat_dir, "dropped_features_log.csv"), index=False)
    
    # Save cleaned merged feature matrix
    cleaned_df = merged_df[["PatientID"] + final_features + ["age", "clinical.T.Stage", "Clinical.N.Stage", "Clinical.M.Stage", "Overall.Stage", "Histology", "gender", "Survival.time", "deadstatus.event"]]
    cleaned_df.to_csv(os.path.join(feat_dir, "cleaned_feature_matrix.csv"), index=False)
    
    return cleaned_df, dropped_df

def calculate_mann_whitney_effect_size(u_stat: float, n1: int, n2: int) -> float:
    """
    Calculates the Rank-Biserial Correlation as the effect size for the Mann-Whitney U test.
    """
    # rank biserial correlation r = 1 - (2 * U) / (n1 * n2)
    # Since U can be defined relative to either sample, we return absolute correlation
    r = 1.0 - (2.0 * u_stat) / (n1 * n2)
    return abs(r)

def run_statistical_analysis(cleaned_df: pd.DataFrame, config: Dict[str, Any]) -> pd.DataFrame:
    """
    Performs univariate association testing between radiomic features and clinical endpoints (Stage 6).
    Applies multiple-testing corrections and reports effect sizes (Rule 18, 19, 20).
    
    Args:
        cleaned_df: Cleaned merged features DataFrame.
        config: Pipeline configuration dictionary.
        
    Returns:
        DataFrame containing statistical outputs for each feature.
    """
    logger.info("Starting Stage 6: Statistical Analysis")
    
    feature_cols = [col for col in cleaned_df.columns if col not in [
        "PatientID", "age", "clinical.T.Stage", "Clinical.N.Stage", 
        "Clinical.M.Stage", "Overall.Stage", "Histology", "gender", 
        "Survival.time", "deadstatus.event"
    ]]
    
    results = []
    
    # Endpoint 1: deadstatus.event (binary)
    dead_group = cleaned_df[cleaned_df["deadstatus.event"] == 1]
    alive_group = cleaned_df[cleaned_df["deadstatus.event"] == 0]
    n_dead = len(dead_group)
    n_alive = len(alive_group)
    
    # Endpoint 2: Overall Stage (ordinal)
    # Map stage to ordinal numbers: I=1, II=2, IIIa=3, IIIb=4
    stage_map = {"I": 1, "II": 2, "IIIa": 3, "IIIb": 4}
    cleaned_df["stage_ordinal"] = cleaned_df["Overall.Stage"].map(stage_map).fillna(0)
    valid_stage_df = cleaned_df[cleaned_df["stage_ordinal"] > 0]
    
    # Endpoint 3: Histology (multi-class)
    # Keep only common classes (squamous cell carcinoma, large cell, adenocarcinoma)
    valid_hist_df = cleaned_df[cleaned_df["Histology"].isin(["squamous cell carcinoma", "large cell", "adenocarcinoma"])]
    
    for feat in feature_cols:
        # --- Survival Status Association (Mann-Whitney U) ---
        u_stat, p_dead = stats.mannwhitneyu(dead_group[feat], alive_group[feat], alternative="two-sided")
        effect_dead = calculate_mann_whitney_effect_size(u_stat, n_dead, n_alive)
        
        # --- Overall Stage Association (Spearman) ---
        rho_stage, p_stage = stats.spearmanr(valid_stage_df[feat], valid_stage_df["stage_ordinal"])
        
        # --- Histology Association (Kruskal-Wallis) ---
        groups = [grp[feat].values for name, grp in valid_hist_df.groupby("Histology")]
        if len(groups) > 1:
            kw_stat, p_hist = stats.kruskal(*groups)
            # Kruskal-Wallis Effect Size: eta^2_H = (H - k + 1) / (N - k)
            k = len(groups)
            N = len(valid_hist_df)
            effect_hist = (kw_stat - k + 1) / (N - k) if N > k else 0.0
        else:
            p_hist = 1.0
            effect_hist = 0.0
            
        results.append({
            "Feature": feat,
            "MW_U_Stat": u_stat,
            "MW_p_value": p_dead,
            "MW_EffectSize_r": effect_dead,
            "Spearman_rho_Stage": rho_stage,
            "Spearman_p_Stage": p_stage,
            "KW_Stat_Hist": kw_stat if len(groups) > 1 else np.nan,
            "KW_p_value": p_hist,
            "KW_EffectSize_eta2": effect_hist
        })
        
    stats_df = pd.DataFrame(results)
    
    # Apply Multiple-Testing Correction (Benjamini-Hochberg FDR) (Rule 18)
    alpha = config["analysis"]["fdr_alpha"]
    
    # Correct survival p-values
    _, fdrs_dead, _, _ = multipletests(stats_df["MW_p_value"], alpha=alpha, method="fdr_bh")
    stats_df["MW_p_FDR"] = fdrs_dead
    
    # Correct stage p-values
    _, fdrs_stage, _, _ = multipletests(stats_df["Spearman_p_Stage"], alpha=alpha, method="fdr_bh")
    stats_df["Spearman_p_FDR"] = fdrs_stage
    
    # Correct histology p-values
    _, fdrs_hist, _, _ = multipletests(stats_df["KW_p_value"], alpha=alpha, method="fdr_bh")
    stats_df["KW_p_FDR"] = fdrs_hist
    
    # Write univariate associations CSV
    table_dir = os.path.join(config["paths"]["output_dir"], "tables")
    stats_df.to_csv(os.path.join(table_dir, "univariate_associations.csv"), index=False)
    
    # Log counts of significant features
    sig_dead = len(stats_df[stats_df["MW_p_FDR"] < alpha])
    sig_stage = len(stats_df[stats_df["Spearman_p_FDR"] < alpha])
    sig_hist = len(stats_df[stats_df["KW_p_FDR"] < alpha])
    
    logger.info(f"Significant associations (FDR < {alpha}): Survival={sig_dead}, Stage={sig_stage}, Histology={sig_hist}")
    return stats_df

def generate_analytical_plots(cleaned_df: pd.DataFrame, stats_df: pd.DataFrame, config: Dict[str, Any]) -> None:
    """
    Generates publication-quality PCA plots, correlation clustermaps, and key feature boxplots (Stage 6).
    """
    logger.info("Generating analytical figures...")
    fig_dir = os.path.join(config["paths"]["output_dir"], "figures")
    os.makedirs(fig_dir, exist_ok=True)
    
    feature_cols = [col for col in cleaned_df.columns if col not in [
        "PatientID", "age", "clinical.T.Stage", "Clinical.N.Stage", 
        "Clinical.M.Stage", "Overall.Stage", "Histology", "gender", 
        "Survival.time", "deadstatus.event", "stage_ordinal"
    ]]
    
    # Standardize features for plotting
    scaler = StandardScaler()
    scaled_feats = scaler.fit_transform(cleaned_df[feature_cols])
    
    # 1. Clustered Heatmap (Spearman Correlation)
    plt.figure(figsize=(12, 10))
    corr_mat = pd.DataFrame(scaled_feats, columns=feature_cols).corr(method="spearman")
    
    g = sns.clustermap(
        corr_mat, 
        cmap="coolwarm", 
        vmin=-1, vmax=1, 
        figsize=(12, 12),
        xticklabels=False, yticklabels=False
    )
    g.fig.suptitle("Spearman Correlation Clustermap of Cleaned Radiomics Features", y=1.02, fontsize=16)
    
    # Include metadata (Rule 17)
    meta_text = f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} | Config v1.0"
    g.fig.text(0.02, 0.02, meta_text, fontsize=9, style="italic")
    
    clustermap_path = os.path.join(fig_dir, "correlation_heatmap.png")
    g.savefig(clustermap_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    # 2. PCA Plot
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled_feats)
    pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"])
    pca_df["Overall.Stage"] = cleaned_df["Overall.Stage"].values
    pca_df["Histology"] = cleaned_df["Histology"].values
    pca_df["SurvivalStatus"] = cleaned_df["deadstatus.event"].map({0: "Alive", 1: "Dead"}).values
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # Plot colored by Stage
    sns.scatterplot(
        data=pca_df[pca_df["Overall.Stage"] != "Unknown"],
        x="PC1", y="PC2",
        hue="Overall.Stage",
        palette="viridis",
        alpha=0.8,
        ax=axes[0]
    )
    axes[0].set_title("PCA colored by Stage", fontsize=12)
    axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)", fontsize=10)
    axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)", fontsize=10)
    
    # Plot colored by Histology
    sns.scatterplot(
        data=pca_df[pca_df["Histology"].isin(["squamous cell carcinoma", "large cell", "adenocarcinoma"])],
        x="PC1", y="PC2",
        hue="Histology",
        palette="Set1",
        alpha=0.8,
        ax=axes[1]
    )
    axes[1].set_title("PCA colored by Histology", fontsize=12)
    axes[1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)", fontsize=10)
    axes[1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)", fontsize=10)
    
    plt.suptitle("Principal Component Analysis (PCA) Projection of Cleaned Radiomics Features", fontsize=14)
    fig.text(0.02, 0.02, meta_text, fontsize=9, style="italic")
    
    pca_path = os.path.join(fig_dir, "pca_plots.png")
    plt.savefig(pca_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    # 3. Box Plots of Top 3 Features associated with Survival Status
    top_survival_feats = stats_df.sort_values(by="MW_p_FDR").head(3)["Feature"].tolist()
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for i, feat in enumerate(top_survival_feats):
        sns.boxplot(
            data=cleaned_df,
            x="deadstatus.event",
            y=feat,
            hue="deadstatus.event",
            palette="Set2",
            ax=axes[i],
            legend=False
        )
        axes[i].set_xticklabels(["Alive (0)", "Dead (1)"])
        axes[i].set_xlabel("Survival Status", fontsize=10)
        axes[i].set_ylabel(feat, fontsize=10)
        
        # Get statistics for title
        row = stats_df[stats_df["Feature"] == feat].iloc[0]
        axes[i].set_title(
            f"{feat.split('_')[-1]}\nFDR p: {row['MW_p_FDR']:.3e} | r: {row['MW_EffectSize_r']:.3f}", 
            fontsize=11
        )
        
    plt.suptitle("Top 3 Radiomic Features Associated with Overall Survival Status", fontsize=14)
    fig.text(0.02, 0.02, meta_text, fontsize=9, style="italic")
    
    box_path = os.path.join(fig_dir, "feature_boxplots_survival.png")
    plt.savefig(box_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"Analytical figures written to {fig_dir}")

def run_analysis_pipeline(config_path: str = "src/config.yaml") -> None:
    """
    Executes Stage 5 (Data Cleaning) and Stage 6 (Statistical Analysis) end-to-end.
    """
    config = load_config(config_path)
    output_dir = config["paths"]["output_dir"]
    raw_features_path = os.path.join(output_dir, "features", "raw_features_all_patients.csv")
    clinical_csv_path = config["paths"]["clinical_csv"]
    
    # Set seed for reproducible PCA or scaling (Rule 16)
    np.random.seed(config["analysis"]["random_seed"])
    
    # Stage 5
    cleaned_df, _ = clean_and_prepare_features(raw_features_path, clinical_csv_path, config)
    
    # Stage 6
    stats_df = run_statistical_analysis(cleaned_df, config)
    generate_analytical_plots(cleaned_df, stats_df, config)
