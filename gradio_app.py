import os
import yaml
import logging
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import joblib
import gradio as gr
from lifelines import KaplanMeierFitter

from src.model_utils import (
    format_ordinal,
    format_risk_markdown,
    get_friendly_feature_name,
    generate_roi_image,
    process_patient_directory
)

# Configure logger
logger = logging.getLogger("radiomics_gradio")
logging.basicConfig(level=logging.INFO)

# Check if running on Hugging Face Spaces (Task 5)
IS_HF_SPACE = (os.environ.get("SYSTEM_ENV") == "HF_SPACE") or (os.environ.get("HF_SPACE") == "true")

# =====================================================================
# 1. LOAD SERIALIZED MODELS & SCALERS (STARTUP)
# =====================================================================

COHORT_PATH = "outputs/features/cleaned_feature_matrix.csv"
MEANINGS_PATH = "outputs/tables/feature_meanings.csv"
SCALER_PATH = "outputs/features/scaler.joblib"
MODEL_RAD_PATH = "outputs/features/model_radiomics.joblib"
MODEL_COMB_PATH = "outputs/features/model_combined.joblib"

if not all(os.path.exists(p) for p in [COHORT_PATH, MEANINGS_PATH, SCALER_PATH, MODEL_RAD_PATH, MODEL_COMB_PATH]):
    raise FileNotFoundError(
        "Required pipeline output files or model checkpoints not found! "
        "Please run the survival pipeline stage first to generate them."
    )

df_cohort = pd.read_csv(COHORT_PATH)
df_cohort = df_cohort.dropna(subset=["Survival.time", "deadstatus.event"])

# Create binary stage and gender indicators
df_cohort["stage_binary"] = df_cohort["Overall.Stage"].apply(lambda s: 1 if s in ["IIIa", "IIIb"] else 0)
df_cohort["gender_binary"] = df_cohort["gender"].map({"male": 1, "female": 0}).fillna(0)
clinical_vars = ["age", "gender_binary", "stage_binary"]

# Load check-pointed objects
scaler = joblib.load(SCALER_PATH)
cph_radiomics = joblib.load(MODEL_RAD_PATH)
cph_combined = joblib.load(MODEL_COMB_PATH)

selected_features = list(cph_radiomics.params_.index)

# Compute PrognosticScore (partial hazards) for the cohort
df_cohort["PrognosticScore"] = cph_radiomics.predict_partial_hazard(scaler.transform(df_cohort[selected_features])).values
median_score = df_cohort["PrognosticScore"].median()
df_cohort["RiskGroup"] = df_cohort["PrognosticScore"].apply(lambda s: "High Risk" if s > median_score else "Low Risk")

# Cache some feature range values for simulations
df_meanings = pd.read_csv(MEANINGS_PATH)
feature_ranges = {}
for feat in selected_features:
    friendly_name = get_friendly_feature_name(feat)
    meaning_row = df_meanings[df_meanings["Feature Name"] == friendly_name]
    feature_ranges[feat] = {
        "min": float(df_cohort[feat].min()),
        "max": float(df_cohort[feat].max()),
        "mean": float(df_cohort[feat].mean()),
        "std": float(df_cohort[feat].std(ddof=0)),
        "def": meaning_row["IBSI Definition"].values[0] if len(meaning_row) > 0 else "Custom feature description.",
        "bio": meaning_row["Biological Interpretation"].values[0] if len(meaning_row) > 0 else "Reflects spatial density variation."
    }

# Sorted patient IDs for dropdown select
patient_ids = sorted(df_cohort["PatientID"].tolist())

# =====================================================================
# 2. HELPER FUNCTIONS
# =====================================================================

# (Helper functions generated programmatically inside src/model_utils.py)

# =====================================================================
# 4. PLOTTING & DIAGNOSTIC INTERFACE CALCULATIONS
# =====================================================================

def calculate_score_and_plot(features_dict, age, gender, stage, patient_id="LUNG1-003", custom_ct_dir=None, custom_seg_path=None):
    # Scale features
    input_vals = [features_dict[feat] for feat in selected_features]
    scaled_vals = scaler.transform(pd.DataFrame([input_vals], columns=selected_features))[0]
    
    # Compute PrognosticScore (relative hazard)
    coefs = cph_radiomics.params_.values
    score = np.exp(np.dot(scaled_vals, coefs))
    
    # Determine risk category
    risk_cat = "High Risk" if score > median_score else "Low Risk"
    
    # Calculate hazard deviation from cohort median
    hazard_ratio_median = score / median_score
    percent_diff = (hazard_ratio_median - 1) * 100
    if percent_diff >= 0:
        hazard_str = f"+{percent_diff:.1f}% above cohort median"
    else:
        hazard_str = f"{percent_diff:.1f}% below cohort median"
        
    # Calculate percentile rank
    percentile_rank = (df_cohort["PrognosticScore"] <= score).mean() * 100
    percentile_str = format_ordinal(percentile_rank)
    
    # Format HTML output
    risk_html = format_risk_markdown(risk_cat, hazard_str, percentile_str)
    
    # Make prediction DataFrame
    patient_df = pd.DataFrame([{
        "PrognosticScore": score,
        "age": age,
        "gender_binary": 1 if gender == "Male" else 0,
        "stage_binary": 1 if stage in ["IIIa", "IIIb"] else 0
    }])
    
    # Get predicted survival curve
    surv_fn = cph_combined.predict_survival_function(patient_df)
    times = surv_fn.index
    survival_probs = surv_fn.values.ravel()
    
    # Plot Kaplan-Meier curves with patient curve overlay
    fig, ax = plt.subplots(figsize=(8, 4.5))
    
    # High-Risk KM
    kmf_high = KaplanMeierFitter()
    high_mask = (df_cohort["RiskGroup"] == "High Risk")
    kmf_high.fit(
        df_cohort.loc[high_mask, "Survival.time"], 
        df_cohort.loc[high_mask, "deadstatus.event"], 
        label="High Risk Reference Cohort"
    )
    kmf_high.plot_survival_function(ax=ax, color="darkred", linestyle="--", ci_show=False, alpha=0.5)
    
    # Low-Risk KM
    kmf_low = KaplanMeierFitter()
    low_mask = (df_cohort["RiskGroup"] == "Low Risk")
    kmf_low.fit(
        df_cohort.loc[low_mask, "Survival.time"], 
        df_cohort.loc[low_mask, "deadstatus.event"], 
        label="Low Risk Reference Cohort"
    )
    kmf_low.plot_survival_function(ax=ax, color="darkgreen", linestyle="--", ci_show=False, alpha=0.5)
    
    # Overlay patient predicted survival curve
    ax.plot(times, survival_probs, color="#38bdf8", linewidth=3.0, label="Patient Predicted Trajectory")
    
    ax.set_title("Prognostic Survival Curve Comparison", fontsize=12, fontweight="bold", color="#f1f5f9")
    ax.set_xlabel("Survival Time (Days)", fontsize=10, color="#f1f5f9")
    ax.set_ylabel("Survival Probability", fontsize=10, color="#f1f5f9")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, df_cohort["Survival.time"].max())
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend(loc="lower left", fontsize=9)
    
    # Set background color to match Gradio dark theme
    fig.patch.set_facecolor("#1e293b")
    ax.set_facecolor("#0f172a")
    ax.tick_params(colors="#f1f5f9")
    for spine in ax.spines.values():
        spine.set_color("#334155")
        
    plt.tight_layout()
    
    # Save curves plot to PNG
    cache_dir = "outputs/figures/cache"
    os.makedirs(cache_dir, exist_ok=True)
    curve_path = os.path.join(cache_dir, f"curve_{patient_id}.png")
    fig.savefig(curve_path, dpi=150, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close(fig)
    
    # Resolve GTV Overlay Visual path
    roi_path = generate_roi_image(patient_id, ct_dir=custom_ct_dir, seg_path=custom_seg_path)
        
    # Format features data table (using clean display names)
    table_rows = []
    for idx, feat in enumerate(selected_features):
        raw = features_dict[feat]
        scaled = scaled_vals[idx]
        definition = feature_ranges[feat]["def"]
        interpretation = feature_ranges[feat]["bio"]
        
        table_rows.append({
            "Feature Name": get_friendly_feature_name(feat),
            "Raw Value": f"{raw:.4f}",
            "Z-Score": f"{scaled:+.3f}",
            "IBSI Definition": definition,
            "Biological Meaning": interpretation
        })
    table_df = pd.DataFrame(table_rows)
    
    return score, risk_html, roi_path, curve_path, table_df

# =====================================================================
# 5. GRADIO APP INTERFACE LAYOUT & STYLING
# =====================================================================

custom_css = """
body { background-color: #0f172a; color: #f1f5f9; }
.gradio-container { max-width: 1200px !important; }
.card { background-color: #1e293b; border: 1px solid #334155; border-radius: 8px; padding: 16px; margin-bottom: 16px; }
"""

with gr.Blocks(title="NSCLC Quantitative Radiomics Signature Demonstrator", css=custom_css) as demo:
    
    gr.HTML("""
    <div style='text-align: center; margin-bottom: 24px;'>
        <h1 style='color: #38bdf8; font-weight: bold; margin-bottom: 4px;'>NSCLC Quantitative Radiomics Signature Demonstrator</h1>
        <p style='color: #94a3b8; font-size: 14px;'>Interactive research demonstrator mapping CT Gross Tumor Volume (GTV) texture phenotypes to Overall Survival outcomes.</p>
    </div>
    """)
    
    with gr.Row():
        # Left Panel - Clinical Context and Diagnostic Outputs
        with gr.Column(scale=5):
            with gr.Group():
                gr.HTML("<h3 style='color:#38bdf8; margin: 4px 0;'>1. Patient Details & Inputs</h3>")
                with gr.Row():
                    input_age = gr.Textbox(label="Age", interactive=False)
                    input_gender = gr.Textbox(label="Gender", interactive=False)
                    input_stage = gr.Textbox(label="Clinical Stage", interactive=False)
                    
            with gr.Group():
                gr.HTML("<h3 style='color:#38bdf8; margin: 4px 0;'>2. Diagnostic Outputs</h3>")
                with gr.Row():
                    out_score = gr.Number(label="Prognostic Score (Relative Hazard)")
                out_risk_html = gr.HTML()
                
            # Model Validation Side Card
            gr.HTML("""
            <div style='background-color: #1e293b; border: 1px solid #334155; border-radius: 8px; padding: 12px; margin-top: 12px; margin-bottom: 12px;'>
                <h3 style='color: #38bdf8; margin: 0 0 8px 0; font-size: 13px; font-weight: bold;'>Model Validation (Lung1 Cohort)</h3>
                <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 8px; font-size: 12px; color: #e2e8f0;'>
                    <div>C-Index: <strong style='color: #38bdf8;'>0.640</strong></div>
                    <div>1-Year AUC: <strong style='color: #38bdf8;'>0.700</strong></div>
                    <div>3-Year AUC: <strong style='color: #38bdf8;'>0.730</strong></div>
                    <div>5-Year AUC: <strong style='color: #38bdf8;'>0.714</strong></div>
                </div>
                <div style='margin-top: 8px; font-size: 11px; color: #94a3b8; border-top: 1px solid #334155; padding-top: 6px;'>
                    Validated on <strong>421</strong> patients via 1000 bootstrap iterations.
                </div>
            </div>
            """)
            
            # CT + Segmentation overlay visualizer
            out_roi_image = gr.Image(label="GTV Segmentation Overlay (Axial CT Zoom)", type="filepath")
            
            # Prognostic survival curves
            out_plot = gr.Image(label="Prognostic Survival Curves", type="filepath")
            
        # Right Panel - Features Table & Interactive Inputs/File Uploads
        with gr.Column(scale=7):
            with gr.Tabs():
                
                # Tab 1: Load Patient Scan
                with gr.TabItem("Load Patient Scan"):
                    gr.HTML("<p style='color:#94a3b8; margin-bottom:8px;'>Select a patient ID from the dropdown, or enter a custom patient ID or full local directory path at the side to extract GTV features.</p>")
                    with gr.Row():
                        select_pid = gr.Dropdown(choices=patient_ids, value="LUNG1-003", label="Select Patient ID")
                        patient_input = gr.Textbox(
                            label="Or Enter Custom Patient ID / Folder Path", 
                            placeholder="e.g. LUNG1-105 or full local directory path...",
                            value="",
                            visible=not IS_HF_SPACE
                        )
                    btn_load = gr.Button("Load Patient Data & Scan", variant="primary")
                    
                # Tab 2: Interactive Simulation Sliders
                with gr.TabItem("Interactive Feature Sliders"):
                    gr.HTML("<p style='color:#94a3b8; margin-bottom:8px;'>Manually adjust selected radiomic features to simulate changes in prognostic outcome.</p>")
                    
                    slider_widgets = {}
                    # Group features for cleaner slider UI
                    with gr.Accordion("Shape & First-Order Features", open=True):
                        for feat in selected_features:
                            if "shape" in feat or "firstorder" in feat:
                                r = feature_ranges[feat]
                                slider_widgets[feat] = gr.Slider(
                                    minimum=r["min"], maximum=r["max"], value=r["mean"], 
                                    label=f"{get_friendly_feature_name(feat)} (Range: {r['min']:.2f} to {r['max']:.2f})"
                                )
                                
                    with gr.Accordion("Texture & Wavelet Features", open=False):
                        for feat in selected_features:
                            if "shape" not in feat and "firstorder" not in feat:
                                r = feature_ranges[feat]
                                slider_widgets[feat] = gr.Slider(
                                    minimum=r["min"], maximum=r["max"], value=r["mean"], 
                                    label=f"{get_friendly_feature_name(feat)} (Range: {r['min']:.2f} to {r['max']:.2f})"
                                )
                    
                    btn_run_sim = gr.Button("Re-Run Simulation Prognosis", variant="secondary")
                    
            # Output Data Table of 19 Selected Features (Dynamic)
            gr.HTML("<h3 style='color:#38bdf8; margin: 12px 0 4px 0;'>3. Feature Interpretations & Raw Values</h3>")
            out_table = gr.Dataframe(
                headers=["Feature Name", "Raw Value", "Z-Score", "IBSI Definition", "Biological Meaning"],
                datatype=["str", "str", "str", "str", "str"],
                wrap=True
            )

    # =====================================================================
    # 6. CALLBACK HANDLERS
    # =====================================================================

    def load_patient_handler(select_pid, patient_input):
        """Unified loader for patient IDs and folder directories."""
        input_str = patient_input.strip() if patient_input else ""
        if not input_str:
            input_str = select_pid if select_pid else ""
            
        input_str = input_str.strip()
        if not input_str:
            return [0.0, "<div style='color:red;'>Error: Please select a Patient ID or enter a Folder Path.</div>", None, None, pd.DataFrame(), "", "", ""] + [gr.update() for _ in selected_features]
            
        ct_dir = None
        seg_path = None
        pid = input_str
        
        # 1. Resolve local path if input is directory
        if os.path.isdir(input_str):
            if IS_HF_SPACE:
                return [0.0, "<div style='color:red;'>Error: Local directory path loading is disabled on public spaces.</div>", None, None, pd.DataFrame(), "", "", ""] + [gr.update() for _ in selected_features]
            from src.data_ingestion import find_patient_series
            discovered_ct, discovered_seg, _ = find_patient_series(input_str)
            if discovered_ct is None or discovered_seg is None:
                return [0.0, f"<div style='color:red;'>Error: Could not locate CT slices and GTV SEG file under: {input_str}</div>", None, None, pd.DataFrame(), "", "", ""] + [gr.update() for _ in selected_features]
                
            ct_dir = discovered_ct
            seg_path = discovered_seg
            pid = os.path.basename(os.path.normpath(input_str))
            
        # 2. Check clinical database
        clinical_rows = df_cohort[df_cohort["PatientID"] == pid]
        if len(clinical_rows) > 0:
            patient_row = clinical_rows.iloc[0]
            age = float(patient_row["age"])
            gender = "Male" if patient_row["gender"] == "male" else "Female"
            stage = str(patient_row["Overall.Stage"])
            
            # Pre-extracted features
            feat_dict = {}
            for feat in selected_features:
                feat_dict[feat] = float(patient_row[feat])
                
            status_msg = f"Success: Loaded data for '{pid}' from database."
        else:
            # If not in database, attempt local path scan in dataset root
            if ct_dir is None or seg_path is None:
                dataset_root = "dataset/NSCLC-Radiomics"
                patient_dir = os.path.join(dataset_root, pid)
                if os.path.isdir(patient_dir):
                    from src.data_ingestion import find_patient_series
                    discovered_ct, discovered_seg, _ = find_patient_series(patient_dir)
                    if discovered_ct is not None and discovered_seg is not None:
                        ct_dir = discovered_ct
                        seg_path = discovered_seg
                        
            if ct_dir is None or seg_path is None:
                return [0.0, f"<div style='color:red;'>Error: Patient '{pid}' not in database and no local DICOM files found.</div>", None, None, pd.DataFrame(), "", "", ""] + [gr.update() for _ in selected_features]
                
            # Extract features on-the-fly
            status_msg, feat_dict = process_patient_directory(ct_dir, seg_path, selected_features)
            if "Success" not in status_msg:
                return [0.0, f"<div style='color:red;'>{status_msg}</div>", None, None, pd.DataFrame(), "", "", ""] + [gr.update() for _ in selected_features]
                
            # Fallback metadata for unrecorded patients
            age = 68.0
            gender = "Male"
            stage = "IIIb"
            status_msg = f"Success: Loaded local scan '{pid}' and extracted GTV features. (Clinical details defaulted)"

        # Calculate prognostic scores, KM curves, and slice crop overlays
        score, risk_html, roi_path, curve_path, table = calculate_score_and_plot(
            feat_dict, age, gender, stage, patient_id=pid, custom_ct_dir=ct_dir, custom_seg_path=seg_path
        )
        
        # Display the load status message above outputs
        styled_status = f"<div style='background-color:#1e293b; padding:8px; border-radius:6px; border:1px solid #334155; margin-bottom:10px; font-size:12px; color:#38bdf8;'>{status_msg}</div>"
        risk_html = styled_status + risk_html
        
        slider_updates = [gr.update(value=feat_dict[feat]) for feat in selected_features]
        
        return [score, risk_html, roi_path, curve_path, table, f"{age:.0f}", gender, stage] + slider_updates

    # Interactive Sliders simulation callback
    def run_slider_simulation(age, gender, stage, *slider_vals):
        feat_dict = {}
        for idx, feat in enumerate(selected_features):
            feat_dict[feat] = slider_vals[idx]
            
        try:
            age_val = float(age)
        except ValueError:
            age_val = 68.0
            
        score, risk_html, roi_path, curve_path, table = calculate_score_and_plot(
            feat_dict, age_val, gender, stage, patient_id="simulated"
        )
        return score, risk_html, curve_path, table

    # Define inputs/outputs lists dynamically
    slider_list = [slider_widgets[feat] for feat in selected_features]

    # Wire Load button
    btn_load.click(
        fn=load_patient_handler,
        inputs=[select_pid, patient_input],
        outputs=[out_score, out_risk_html, out_roi_image, out_plot, out_table, input_age, input_gender, input_stage] + slider_list
    )

    # Wire Simulation button
    btn_run_sim.click(
        fn=run_slider_simulation,
        inputs=[input_age, input_gender, input_stage] + slider_list,
        outputs=[out_score, out_risk_html, out_plot, out_table]
    )

    # Trigger loading patient LUNG1-003 on launch
    demo.load(
        fn=load_patient_handler,
        inputs=[select_pid, patient_input],
        outputs=[out_score, out_risk_html, out_roi_image, out_plot, out_table, input_age, input_gender, input_stage] + slider_list
    )

# Launch
if __name__ == "__main__":
    username = os.environ.get("GRADIO_AUTH_USER")
    password = os.environ.get("GRADIO_AUTH_PASS")
    auth_creds = (username, password) if username and password else None
    
    if auth_creds:
        logger.info("Launching with basic authentication enabled.")
        demo.launch(server_name="0.0.0.0", share=False, auth=auth_creds)
    else:
        demo.launch(server_name="0.0.0.0", share=False)

