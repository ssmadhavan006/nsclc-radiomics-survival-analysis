import os
import glob
import pandas as pd
import pydicom
import logging
from datetime import datetime
from typing import Dict, Tuple, List, Optional

logger = logging.getLogger("radiomics_pipeline")

def load_clinical_data(csv_path: str) -> pd.DataFrame:
    """
    Loads clinical metadata from CSV and performs basic validation.
    
    Args:
        csv_path: Path to the clinical CSV file.
        
    Returns:
        A validated pandas DataFrame containing clinical data.
        
    Raises:
        ValueError: If essential columns are missing.
    """
    logger.info(f"Loading clinical data from {csv_path}")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Clinical CSV not found at: {csv_path}")
        
    df = pd.read_csv(csv_path)
    
    # Required columns for our pipeline
    required_cols = [
        "PatientID", "age", "clinical.T.Stage", "Clinical.N.Stage", 
        "Clinical.M.Stage", "Overall.Stage", "Histology", "gender", 
        "Survival.time", "deadstatus.event"
    ]
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Clinical CSV is missing required columns: {missing_cols}")
        
    # Standardize types and fill simple categories
    df["PatientID"] = df["PatientID"].astype(str).str.strip()
    df["Histology"] = df["Histology"].fillna("Unknown").astype(str).str.strip()
    df["Overall.Stage"] = df["Overall.Stage"].fillna("Unknown").astype(str).str.strip()
    
    logger.info(f"Successfully loaded clinical data for {len(df)} patients.")
    return df

def detect_dicom_modality(filepath: str) -> str:
    """
    Reads the metadata of a DICOM file to determine if it is SEG, RTSTRUCT, or something else.
    
    Args:
        filepath: Path to the DICOM file.
        
    Returns:
        String indicator of modality: 'SEG', 'RTSTRUCT', 'CT', or 'UNKNOWN'
    """
    try:
        header = pydicom.dcmread(filepath, stop_before_pixels=True)
        sop_class_uid = getattr(header, "SOPClassUID", "")
        modality = getattr(header, "Modality", "")
        
        # Standard SOP Class UIDs
        # SEG: 1.2.840.10008.5.1.4.1.1.66.4
        # RTSTRUCT: 1.2.840.10008.5.1.4.1.1.481.3
        # CT: 1.2.840.10008.5.1.4.1.1.2
        
        if sop_class_uid == "1.2.840.10008.5.1.4.1.1.66.4" or modality == "SEG":
            return "SEG"
        elif sop_class_uid == "1.2.840.10008.5.1.4.1.1.481.3" or modality == "RTSTRUCT":
            return "RTSTRUCT"
        elif sop_class_uid == "1.2.840.10008.5.1.4.1.1.2" or modality == "CT":
            return "CT"
        else:
            return str(modality) if modality else "UNKNOWN"
            
    except Exception as e:
        logger.warning(f"Failed to read DICOM header for {filepath}: {str(e)}")
        return "ERROR"

def find_patient_series(patient_dir: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Walks a patient directory to identify the CT folder, the Segmentation folder, and Seg file.
    
    Args:
        patient_dir: Absolute path to the patient's directory.
        
    Returns:
        Tuple: (ct_folder_path, seg_file_path, seg_modality)
    """
    # Enumerate the single study directory
    subdirs = [os.path.join(patient_dir, d) for d in os.listdir(patient_dir) 
               if os.path.isdir(os.path.join(patient_dir, d))]
    if not subdirs:
        return None, None, None
        
    study_dir = subdirs[0]
    series_dirs = [os.path.join(study_dir, s) for s in os.listdir(study_dir) 
                   if os.path.isdir(os.path.join(study_dir, s))]
                   
    ct_dir = None
    seg_file = None
    seg_modality = None
    
    for s_dir in series_dirs:
        basename = os.path.basename(s_dir)
        dcm_files = glob.glob(os.path.join(s_dir, "*.dcm"))
        
        if not dcm_files:
            continue
            
        # Segmentation check: starts with 300 or name contains Segmentation
        if basename.startswith("300.") or "segmentation" in basename.lower():
            if len(dcm_files) == 1:
                test_file = dcm_files[0]
                modality = detect_dicom_modality(test_file)
                if modality in ["SEG", "RTSTRUCT"]:
                    seg_file = test_file
                    seg_modality = modality
                else:
                    logger.debug(f"File {test_file} has modality {modality}, skipped as segmentation candidate.")
            else:
                logger.warning(f"Segmentation folder {s_dir} contains multiple files: {len(dcm_files)}")
        else:
            # CT check: should have multiple dcm slices
            if len(dcm_files) > 1:
                # Confirm it's actually CT
                modality = detect_dicom_modality(dcm_files[0])
                if modality == "CT":
                    ct_dir = s_dir
                else:
                    logger.debug(f"Directory {s_dir} contains multi-slice DICOMs but modality is {modality}")
                    
    return ct_dir, seg_file, seg_modality

def build_data_manifest(dataset_root: str, output_dir: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Scans the dataset folder and clinical CSV to construct a data manifest.
    Tracks failures and writes them to failed_cases.csv per Rule 7.
    
    Args:
        dataset_root: Path containing patient folders.
        output_dir: Main output directory for saving logs.
        
    Returns:
        Tuple: (manifest_df, failed_cases_df)
    """
    logger.info(f"Scanning patient folders under {dataset_root}")
    patient_ids = [d for d in os.listdir(dataset_root) 
                   if os.path.isdir(os.path.join(dataset_root, d)) and d.startswith("LUNG1-")]
    
    manifest_records: List[Dict[str, Any]] = []
    failed_records: List[Dict[str, Any]] = []
    timestamp_str = datetime.now().isoformat()
    
    for pid in sorted(patient_ids):
        p_dir = os.path.join(dataset_root, pid)
        ct_dir, seg_file, seg_modality = find_patient_series(p_dir)
        
        failure_reasons = []
        if not ct_dir:
            failure_reasons.append("Missing CT series")
        if not seg_file:
            failure_reasons.append("Missing Segmentation series (SEG/RTSTRUCT)")
            
        if failure_reasons:
            reason = "; ".join(failure_reasons)
            logger.warning(f"Patient {pid} failed ingestion checks: {reason}")
            failed_records.append({
                "PatientID": pid,
                "Stage": "Ingestion",
                "Reason": reason,
                "Timestamp": timestamp_str
            })
            manifest_records.append({
                "PatientID": pid,
                "CTPath": ct_dir if ct_dir else "",
                "SegPath": seg_file if seg_file else "",
                "SegModality": seg_modality if seg_modality else "",
                "Status": "invalid",
                "FailureReason": reason
            })
        else:
            manifest_records.append({
                "PatientID": pid,
                "CTPath": ct_dir,
                "SegPath": seg_file,
                "SegModality": seg_modality,
                "Status": "valid",
                "FailureReason": ""
            })
            
    manifest_df = pd.DataFrame(manifest_records)
    failed_df = pd.DataFrame(failed_records)
    
    # Save manifest
    logs_dir = os.path.join(output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    
    manifest_path = os.path.join(logs_dir, "data_manifest.csv")
    manifest_df.to_csv(manifest_path, index=False)
    logger.info(f"Manifest written to {manifest_path} ({len(manifest_df)} patients, {len(manifest_df[manifest_df['Status'] == 'valid'])} valid)")
    
    # Save failed cases to failed_cases.csv in outputs/logs/
    failed_cases_path = os.path.join(logs_dir, "failed_cases.csv")
    failed_df.to_csv(failed_cases_path, index=False)
    logger.info(f"Failed cases logged to {failed_cases_path} ({len(failed_df)} cases)")
    
    return manifest_df, failed_df
