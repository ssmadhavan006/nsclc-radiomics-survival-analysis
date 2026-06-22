import os
import zipfile
import tempfile
import shutil
import yaml
import logging
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import gradio as gr
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter, CoxPHFitter

# Configure logger
logger = logging.getLogger("radiomics_gradio")
logging.basicConfig(level=logging.INFO)

# =====================================================================
# 1. MODEL FITTING & SCALER SETUP ON PIPELINE OUTPUTS (STARTUP)
# =====================================================================

# Load the cohort features and outcomes
COHORT_PATH = "outputs/features/cleaned_feature_matrix.csv"
MEANINGS_PATH = "outputs/tables/feature_meanings.csv"

if not os.path.exists(COHORT_PATH) or not os.path.exists(MEANINGS_PATH):
    raise FileNotFoundError(
        "Required pipeline output files not found! Please run the survival pipeline stage first to generate: "
        f"{COHORT_PATH} and {MEANINGS_PATH}"
    )

df_cohort = pd.read_csv(COHORT_PATH)
df_cohort = df_cohort.dropna(subset=["Survival.time", "deadstatus.event"])

# Create binary stage and gender indicators
df_cohort["stage_binary"] = df_cohort["Overall.Stage"].apply(lambda s: 1 if s in ["IIIa", "IIIb"] else 0)
df_cohort["gender_binary"] = df_cohort["gender"].map({"male": 1, "female": 0}).fillna(0)
clinical_vars = ["age", "gender_binary", "stage_binary"]

# Load selected 19 features
df_meanings = pd.read_csv(MEANINGS_PATH)
selected_features = df_meanings["Feature Name"].tolist()

# Fit StandardScaler on cohort features
scaler = StandardScaler()
scaler.fit(df_cohort[selected_features])

# Scale selected features to fit radiomics signature
scaled_features = pd.DataFrame(scaler.transform(df_cohort[selected_features]), columns=selected_features)
scaled_features["Survival.time"] = df_cohort["Survival.time"].values
scaled_features["deadstatus.event"] = df_cohort["deadstatus.event"].values

# Fit CPH Radiomics Signature model to predict PrognosticScore
cph_radiomics = CoxPHFitter(penalizer=0.05, l1_ratio=1.0)
cph_radiomics.fit(
    scaled_features[selected_features + ["Survival.time", "deadstatus.event"]], 
    "Survival.time", 
    "deadstatus.event"
)

# Compute PrognosticScore (partial hazards) for the cohort
df_cohort["PrognosticScore"] = cph_radiomics.predict_partial_hazard(scaled_features[selected_features]).values
median_score = df_cohort["PrognosticScore"].median()
df_cohort["RiskGroup"] = df_cohort["PrognosticScore"].apply(lambda s: "High Risk" if s > median_score else "Low Risk")

# Fit combined multivariate Cox model
cph_combined = CoxPHFitter()
cph_combined.fit(
    df_cohort[["PrognosticScore"] + clinical_vars + ["Survival.time", "deadstatus.event"]], 
    "Survival.time", 
    "deadstatus.event"
)

# Cache some feature range values for simulations
feature_ranges = {}
for feat in selected_features:
    feature_ranges[feat] = {
        "min": float(df_cohort[feat].min()),
        "max": float(df_cohort[feat].max()),
        "mean": float(df_cohort[feat].mean()),
        "std": float(df_cohort[feat].std(ddof=0)),
        "def": df_meanings.loc[df_meanings["Feature Name"] == feat, "IBSI Definition"].values[0],
        "bio": df_meanings.loc[df_meanings["Feature Name"] == feat, "Biological Interpretation"].values[0]
    }

# Sorted patient IDs for dropdown select
patient_ids = sorted(df_cohort["PatientID"].tolist())

# =====================================================================
# 2. HELPER FUNCTIONS
# =====================================================================

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
        <strong style='color: {color}; font-size: 22px;'>{risk_cat}</strong>
    </div>
    <div style='margin-bottom: 12px; border-top: 1px solid #334155; padding-top: 8px;'>
        <span style='color: #94a3b8; font-size: 11px; text-transform: uppercase; font-weight: bold;'>Predicted Hazard Deviation</span><br/>
        <strong style='color: #f1f5f9; font-size: 15px;'>{hazard_str}</strong>
    </div>
    <div style='border-top: 1px solid #334155; padding-top: 8px;'>
        <span style='color: #94a3b8; font-size: 11px; text-transform: uppercase; font-weight: bold;'>Cohort Percentile Rank</span><br/>
        <strong style='color: #38bdf8; font-size: 15px;'>{percentile_str}</strong>
    </div>
</div>
"""

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
            # Fallback to middle slice
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

# =====================================================================
# 3. PIPELINE UPLOADS PROCESSING LOGIC
# =====================================================================

def process_uploaded_dicoms(zip_file, seg_file):
    if zip_file is None or seg_file is None:
        return "Error: Please upload both the CT slices (.zip) and the GTV mask (.dcm)", None
        
    temp_dir = tempfile.mkdtemp()
    try:
        # Extract CT ZIP
        with zipfile.ZipFile(zip_file.name, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)
            
        # Discover CT files (ignore the SEG file if placed inside)
        seg_name = os.path.basename(seg_file.name)
        dcm_files = []
        for root, _, files in os.walk(temp_dir):
            for f in files:
                if f.endswith('.dcm') and f != seg_name:
                    dcm_files.append(os.path.join(root, f))
                    
        if not dcm_files:
            return "Error: No CT slice DICOM files (.dcm) found in the uploaded ZIP.", None
            
        ct_dir = os.path.dirname(dcm_files[0])
        
        # Load master config
        with open("src/config.yaml", "r") as f:
            config = yaml.safe_load(f)
            
        # Run ingestion & preprocessing steps
        from src.preprocessing import preprocess_case
        preprocessed_ct, preprocessed_mask = preprocess_case(ct_dir, seg_file.name, config)
        
        # Check overlaps
        mask_arr = sitk.GetArrayFromImage(preprocessed_mask)
        voxel_count = np.sum(mask_arr)
        if voxel_count < 50:
            return f"Error: Resampled GTV mask only has {voxel_count} voxels. Must contain at least 50.", None
            
        # Run PyRadiomics feature extractor
        from src.feature_extraction import get_radiomics_extractor
        extractor = get_radiomics_extractor(config)
        feature_vector = extractor.execute(preprocessed_ct, preprocessed_mask)
        
        # Filter out 19 features
        extracted_19 = {}
        for feat in selected_features:
            val = feature_vector.get(feat, 0.0)
            if hasattr(val, "item"):
                extracted_19[feat] = float(val.item())
            else:
                extracted_19[feat] = float(val)
                
        return "Success: CT volume and GTV segmentation successfully matched and processed!", extracted_19
        
    except Exception as e:
        return f"Pipeline execution failed: {str(e)}", None
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

# =====================================================================
# 4. PLOTTING & DIAGNOSTIC INTERFACE CALCULATIONS
# =====================================================================

def calculate_score_and_plot(features_dict, age, gender, stage, patient_id="LUNG1-003", custom_roi_path=None):
    # Scale features
    input_vals = [features_dict[feat] for feat in selected_features]
    scaled_vals = scaler.transform([input_vals])[0]
    
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
    
    # Resolve the ROI overlay path
    if custom_roi_path:
        roi_path = custom_roi_path
    else:
        roi_path = generate_roi_image(patient_id)
        
    # Format features data table (using distinct names to prevent name collision)
    table_rows = []
    for idx, feat in enumerate(selected_features):
        raw = features_dict[feat]
        scaled = scaled_vals[idx]
        definition = feature_ranges[feat]["def"]
        interpretation = feature_ranges[feat]["bio"]
        
        # Display name includes wavelet and modality prefix for uniqueness
        parts = feat.split("_")
        if len(parts) >= 3:
            display_name = f"{parts[-3]}_{parts[-2]}_{parts[-1]}"
        elif len(parts) == 2:
            display_name = f"{parts[-2]}_{parts[-1]}"
        else:
            display_name = parts[-1]
            
        display_name = display_name.replace("wavelet-", "")
        
        table_rows.append({
            "Feature Name": display_name,
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
                    input_age = gr.Slider(minimum=30, maximum=95, value=68, step=1, label="Age")
                    input_gender = gr.Dropdown(choices=["Male", "Female"], value="Male", label="Gender")
                    input_stage = gr.Dropdown(choices=["I", "II", "IIIa", "IIIb"], value="IIIb", label="Clinical Stage")
                    
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
                
                # Tab 1: Cohort Select Demo
                with gr.TabItem("Demo Patient Selection"):
                    gr.HTML("<p style='color:#94a3b8; margin-bottom:8px;'>Choose a patient from the validated TCIA Lung1 cohort to inspect their GTV radiomic signature and CT segment overlay.</p>")
                    select_pid = gr.Dropdown(choices=patient_ids, value="LUNG1-003", label="Select Patient ID")
                    btn_load_demo = gr.Button("Load Demo Patient Data", variant="primary")
                    
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
                                    label=f"{feat.split('_')[-1]} (Scale: {r['min']:.1f} to {r['max']:.1f})"
                                )
                                
                    with gr.Accordion("Texture & Wavelet Features", open=False):
                        for feat in selected_features:
                            if "shape" not in feat and "firstorder" not in feat:
                                r = feature_ranges[feat]
                                slider_widgets[feat] = gr.Slider(
                                    minimum=r["min"], maximum=r["max"], value=r["mean"], 
                                    label=f"{feat.split('_')[-2]}_{feat.split('_')[-1]} (Wavelet)"
                                )
                    
                    btn_run_sim = gr.Button("Re-Run Simulation Prognosis", variant="secondary")
                    
                # Tab 3: Custom DICOM Upload
                with gr.TabItem("Custom DICOM Upload"):
                    gr.HTML("<p style='color:#94a3b8; margin-bottom:8px;'>Upload custom pre-treatment CT slices and GTV SEG file to run the coordinates matching and extraction pipeline in real-time.</p>")
                    with gr.Row():
                        upload_zip = gr.File(label="CT Slices (ZIP containing .dcm slices)", file_types=[".zip"])
                        upload_seg = gr.File(label="GTV Mask (DICOM SEG file .dcm)", file_types=[".dcm"])
                    upload_status = gr.Textbox(label="Pipeline Processing Log", placeholder="Awaiting file upload...")
                    btn_run_upload = gr.Button("Run Real-Time Extraction Pipeline", variant="primary")
                    
            # Output Data Table of 19 Selected Features (Dynamic)
            gr.HTML("<h3 style='color:#38bdf8; margin: 12px 0 4px 0;'>3. Feature Interpretations & Raw Values</h3>")
            out_table = gr.Dataframe(
                headers=["Feature Name", "Raw Value", "Z-Score", "IBSI Definition", "Biological Meaning"],
                datatype=["str", "str", "str", "str", "str"],
                wrap=True
            )

    # =====================================================================
    # 6. BUTTON CLICK HANDLERS
    # =====================================================================

    # Dropdown select loader callback
    def load_demo_data(pid):
        patient_row = df_cohort[df_cohort["PatientID"] == pid].iloc[0]
        
        # Extracted features dictionary
        feat_dict = {}
        for feat in selected_features:
            feat_dict[feat] = float(patient_row[feat])
            
        # Clinical parameters
        age = float(patient_row["age"])
        gender = "Male" if patient_row["gender"] == "male" else "Female"
        stage = str(patient_row["Overall.Stage"])
        
        # Calculate C-index and generate curves
        score, risk_html, roi_path, curve_path, table = calculate_score_and_plot(feat_dict, age, gender, stage, patient_id=pid)
        
        # Generate sliders updates
        slider_updates = [gr.update(value=feat_dict[feat]) for feat in selected_features]
        
        return [score, risk_html, roi_path, curve_path, table, age, gender, stage] + slider_updates

    # Interactive Sliders simulation callback (only updates score/curves, doesn't change patient scan)
    def run_slider_simulation(age, gender, stage, *slider_vals):
        feat_dict = {}
        for idx, feat in enumerate(selected_features):
            feat_dict[feat] = slider_vals[idx]
            
        score, risk_html, roi_path, curve_path, table = calculate_score_and_plot(feat_dict, age, gender, stage, patient_id="simulated")
        return score, risk_html, curve_path, table

    # DICOM upload pipeline callback
    def run_upload_pipeline(zip_file, seg_file, age, gender, stage):
        status, feat_dict = process_uploaded_dicoms(zip_file, seg_file)
        if "Success" not in status:
            return status, None, None, None, None, None
            
        patient_id = "custom_upload"
        roi_img_path = generate_roi_image(patient_id, ct_dir=zip_file.name, seg_path=seg_file.name)
        
        score, risk_html, roi_path, curve_path, table = calculate_score_and_plot(
            feat_dict, age, gender, stage, patient_id=patient_id, custom_roi_path=roi_img_path
        )
        
        # Generate sliders updates to match uploaded values
        slider_updates = [gr.update(value=feat_dict[feat]) for feat in selected_features]
        
        return [status, score, risk_html, roi_path, curve_path, table] + slider_updates

    # Define inputs/outputs lists dynamically
    slider_list = [slider_widgets[feat] for feat in selected_features]

    # Wire Load Demo button
    btn_load_demo.click(
        fn=load_demo_data,
        inputs=[select_pid],
        outputs=[out_score, out_risk_html, out_roi_image, out_plot, out_table, input_age, input_gender, input_stage] + slider_list
    )

    # Wire Simulation button
    btn_run_sim.click(
        fn=run_slider_simulation,
        inputs=[input_age, input_gender, input_stage] + slider_list,
        outputs=[out_score, out_risk_html, out_plot, out_table]
    )

    # Wire Custom Upload button
    btn_run_upload.click(
        fn=run_upload_pipeline,
        inputs=[upload_zip, upload_seg, input_age, input_gender, input_stage],
        outputs=[upload_status, out_score, out_risk_html, out_roi_image, out_plot, out_table] + slider_list
    )

    # Trigger loading patient LUNG1-003 on launch
    demo.load(
        fn=load_demo_data,
        inputs=[select_pid],
        outputs=[out_score, out_risk_html, out_roi_image, out_plot, out_table, input_age, input_gender, input_stage] + slider_list
    )

# Launch
if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", share=False)
