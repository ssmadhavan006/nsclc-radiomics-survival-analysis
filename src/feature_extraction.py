import os
import time
import glob
import yaml
import logging
import pandas as pd
from datetime import datetime
from typing import Dict, Any, List, Tuple
from radiomics import featureextractor
import SimpleITK as sitk
from joblib import Parallel, delayed

from src.utils import setup_logger, load_config
from src.preprocessing import preprocess_case

logger = logging.getLogger("radiomics_pipeline")

def get_radiomics_extractor(config: Dict[str, Any]) -> featureextractor.RadiomicsFeatureExtractor:
    """
    Creates and configures a PyRadiomics Feature Extractor programmatically.
    
    Args:
        config: Pipeline configuration dictionary.
        
    Returns:
        Configured RadiomicsFeatureExtractor instance.
    """
    extractor = featureextractor.RadiomicsFeatureExtractor()
    
    # 1. Update settings
    settings = config["radiomics"]["setting"]
    extractor.settings.update(settings)
    
    # 2. Configure Image Types
    extractor.disableAllImageTypes()
    for img_type, params in config["radiomics"]["imageType"].items():
        extractor.enableImageTypeByName(img_type, enabled=True, customArgs=params if params else {})
        
    # 3. Configure Feature Classes
    extractor.disableAllFeatures()
    for feat_class in config["radiomics"]["featureClass"]:
        extractor.enableFeatureClassByName(feat_class)
        
    return extractor

def save_extractor_settings(extractor: featureextractor.RadiomicsFeatureExtractor, output_dir: str) -> None:
    """
    Saves the exact configuration settings used by PyRadiomics for version control (Rule 13).
    
    Args:
        extractor: The configured RadiomicsFeatureExtractor.
        output_dir: Main output directory.
    """
    feat_dir = os.path.join(output_dir, "features")
    os.makedirs(feat_dir, exist_ok=True)
    param_path = os.path.join(feat_dir, "pyradiomics_params.yaml")
    
    config_dict = {
        "settings": extractor.settings,
        "enabledImageTypes": extractor.enabledImagetypes,
        "enabledFeatures": extractor.enabledFeatures
    }
    
    try:
        with open(param_path, "w", encoding="utf-8") as f:
            yaml.safe_dump(config_dict, f, default_flow_style=False)
        logger.info(f"Saved PyRadiomics parameter configuration to {param_path}")
    except Exception as e:
        logger.error(f"Failed to save PyRadiomics parameters: {str(e)}")

def process_single_patient(
    pid: str,
    ct_dir: str,
    seg_path: str,
    config_path: str
) -> Tuple[str, str, float]:
    """
    Processes a single patient: loads images, preprocesses, extracts features, 
    and writes a temporary CSV file. This is executed in parallel processes.
    
    Args:
        pid: Patient ID.
        ct_dir: Directory containing CT slices.
        seg_path: Path to the DICOM SEG file.
        config_path: Path to the config file.
        
    Returns:
        Tuple: (patient_id, status ['success' or 'error'], execution_time)
    """
    start_time = time.time()
    
    # Configure a local logger for the process
    # (Since logging within joblib parallel processes can sometimes be tricky)
    local_logger = logging.getLogger(f"radiomics_pipeline.{pid}")
    
    try:
        config = load_config(config_path)
        output_dir = config["paths"]["output_dir"]
        feat_dir = os.path.join(output_dir, "features")
        logs_dir = os.path.join(output_dir, "logs")
        
        temp_csv_path = os.path.join(feat_dir, f"temp_{pid}.csv")
        temp_fail_path = os.path.join(feat_dir, f"fail_{pid}.csv")
        
        # 1. Check if already successfully processed
        if os.path.exists(temp_csv_path):
            return pid, "success", 0.0
            
        # 2. Run Preprocessing
        preprocessed_ct, preprocessed_mask = preprocess_case(ct_dir, seg_path, config)
        
        # 3. Initialize Extractor & Run Extraction
        extractor = get_radiomics_extractor(config)
        feature_vector = extractor.execute(preprocessed_ct, preprocessed_mask)
        
        # 4. Save to temporary CSV file
        patient_record = {"PatientID": pid}
        for k, v in feature_vector.items():
            if isinstance(v, (int, float)):
                patient_record[k] = float(v)
            elif hasattr(v, "item"):
                patient_record[k] = float(v.item())
            else:
                patient_record[k] = str(v)
                
        pd.DataFrame([patient_record]).to_csv(temp_csv_path, index=False)
        
        # Remove old failure checkpoint if it exists
        if os.path.exists(temp_fail_path):
            os.remove(temp_fail_path)
            
        duration = time.time() - start_time
        return pid, "success", duration
        
    except Exception as e:
        duration = time.time() - start_time
        error_msg = str(e)
        
        # Save failure details
        config = load_config(config_path)
        feat_dir = os.path.join(config["paths"]["output_dir"], "features")
        temp_fail_path = os.path.join(feat_dir, f"fail_{pid}.csv")
        
        fail_record = {
            "PatientID": pid,
            "Stage": "Extraction",
            "Reason": error_msg,
            "Timestamp": datetime.now().isoformat(),
            "Duration": duration
        }
        pd.DataFrame([fail_record]).to_csv(temp_fail_path, index=False)
        
        return pid, "error", duration

def run_extraction_pipeline(config_path: str = "src/config.yaml", limit: int = -1) -> None:
    """
    Runs the feature extraction pipeline in parallel across multiple CPU cores.
    Consolidates temporary outputs into the final raw features matrix.
    
    Args:
        config_path: Path to config.yaml.
        limit: Limit number of cases for testing.
    """
    config = load_config(config_path)
    output_dir = config["paths"]["output_dir"]
    logs_dir = os.path.join(output_dir, "logs")
    feat_dir = os.path.join(output_dir, "features")
    
    manifest_path = os.path.join(logs_dir, "data_manifest.csv")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Data manifest not found at: {manifest_path}. Run Stage 2 first.")
        
    manifest = pd.read_csv(manifest_path)
    valid_cases = manifest[manifest["Status"] == "valid"]
    
    if limit > 0:
        valid_cases = valid_cases.head(limit)
        logger.info(f"Limiting execution to first {limit} patients.")
        
    logger.info(f"Preparing parallel feature extraction for {len(valid_cases)} valid patients.")
    
    # Save the configuration settings for reproducibility
    extractor = get_radiomics_extractor(config)
    save_extractor_settings(extractor, output_dir)
    
    # Build list of tasks
    tasks = []
    skipped_count = 0
    for _, row in valid_cases.iterrows():
        pid = str(row["PatientID"]).strip()
        temp_csv_path = os.path.join(feat_dir, f"temp_{pid}.csv")
        
        if os.path.exists(temp_csv_path):
            skipped_count += 1
            continue
            
        tasks.append((pid, row["CTPath"], row["SegPath"]))
        
    logger.info(f"Skipped {skipped_count} already-processed patients. {len(tasks)} tasks to execute.")
    
    if tasks:
        # Run parallel extraction loop (using all available CPU cores)
        # Using n_jobs=-1 for maximum efficiency
        logger.info(f"Spawning parallel jobs with n_jobs=-1...")
        results = Parallel(n_jobs=-1, verbose=10)(
            delayed(process_single_patient)(pid, ct_dir, seg_path, config_path)
            for pid, ct_dir, seg_path in tasks
        )
        
        # Log summary of results
        succeeded = [pid for pid, status, _ in results if status == "success"]
        failed = [pid for pid, status, _ in results if status == "error"]
        logger.info(f"Batch execution complete. Successes: {len(succeeded)}, Failures: {len(failed)}")
        
    # Consolidate all temporary CSVs into the final raw features matrix
    logger.info("Consolidating temporary patient CSVs...")
    all_temp_files = glob.glob(os.path.join(feat_dir, "temp_LUNG1-*.csv"))
    
    if all_temp_files:
        df_list = []
        for f in sorted(all_temp_files):
            try:
                df_list.append(pd.read_csv(f))
            except Exception as e:
                logger.error(f"Error reading temporary file {f}: {str(e)}")
                
        if df_list:
            raw_features_df = pd.concat(df_list, ignore_index=True)
            raw_features_path = os.path.join(feat_dir, "raw_features_all_patients.csv")
            raw_features_df.to_csv(raw_features_path, index=False)
            logger.info(f"Aggregated raw feature matrix written to {raw_features_path} with shape {raw_features_df.shape}")
            
    # Consolidate all failed cases
    all_fail_files = glob.glob(os.path.join(feat_dir, "fail_LUNG1-*.csv"))
    failed_extractions_path = os.path.join(feat_dir, "failed_extractions.csv")
    global_failed_cases_path = os.path.join(logs_dir, "failed_cases.csv")
    
    if all_fail_files:
        fail_list = []
        for f in sorted(all_fail_files):
            try:
                fail_list.append(pd.read_csv(f))
            except Exception:
                pass
                
        if fail_list:
            fail_df = pd.concat(fail_list, ignore_index=True)
            fail_df.to_csv(failed_extractions_path, index=False)
            logger.info(f"Failed extractions list written to {failed_extractions_path} ({len(fail_df)} failures)")
            
            # Merge with existing failed cases from ingestion
            existing_fails = []
            if os.path.exists(global_failed_cases_path):
                try:
                    existing_fails.append(pd.read_csv(global_failed_cases_path))
                except Exception:
                    pass
            existing_fails.append(fail_df[["PatientID", "Stage", "Reason", "Timestamp"]])
            
            pd.concat(existing_fails, ignore_index=True).drop_duplicates(subset=["PatientID", "Stage"]).to_csv(global_failed_cases_path, index=False)
            logger.info("Updated global failed_cases.csv.")
            
    logger.info("Feature extraction pipeline stage complete.")
