import os
import argparse
import logging
from src.utils import setup_logger, load_config, log_stage
from src.data_ingestion import load_clinical_data, build_data_manifest
from src.feature_extraction import run_extraction_pipeline
from src.analysis import run_analysis_pipeline
from src.survival import run_survival_pipeline

def main() -> None:
    """
    Main entry point for orchestrating the radiomics research pipeline (Rule 15).
    """
    parser = argparse.ArgumentParser(
        description="End-to-End reproducible radiomics research pipeline for NSCLC-Radiomics Lung1."
    )
    parser.add_argument(
        "--config", 
        type=str, 
        default="src/config.yaml", 
        help="Path to the master configuration YAML file."
    )
    parser.add_argument(
        "--stage", 
        type=str, 
        default="all", 
        choices=["ingestion", "extraction", "analysis", "survival", "all"], 
        help="Which pipeline stage to execute."
    )
    parser.add_argument(
        "--limit", 
        type=int, 
        default=-1, 
        help="Limit number of patients to process (for fast testing/debugging)."
    )
    args = parser.parse_args()
    
    # 1. Setup Logging
    logger = setup_logger()
    
    # Load configuration
    try:
        config = load_config(args.config)
    except Exception as e:
        logger.error(f"Failed to load configuration from {args.config}: {str(e)}")
        return
        
    dataset_root = config["paths"]["dataset_root"]
    clinical_csv = config["paths"]["clinical_csv"]
    output_dir = config["paths"]["output_dir"]
    
    logger.info("Initializing NSCLC-Radiomics End-to-End Pipeline")
    logger.info(f"Target Stage: {args.stage.upper()}")
    
    # Run Ingestion Stage
    if args.stage in ["ingestion", "all"]:
        with log_stage(logger, "Data Ingestion"):
            # Load and validate clinical CSV
            clinical_df = load_clinical_data(clinical_csv)
            # Scan patient directories and build manifest
            manifest_df, failed_df = build_data_manifest(dataset_root, output_dir)
            logger.info(f"Ingestion summary: {len(manifest_df)} patients scanned, {len(manifest_df[manifest_df['Status'] == 'valid'])} valid.")
            
    # Run Extraction Stage
    if args.stage in ["extraction", "all"]:
        with log_stage(logger, "Radiomic Feature Extraction"):
            run_extraction_pipeline(args.config, limit=args.limit)
            
    # Run Analysis Stage (Stage 5 + 6)
    if args.stage in ["analysis", "all"]:
        with log_stage(logger, "Data Cleaning & Statistical Analysis"):
            run_analysis_pipeline(args.config)
            
    # Run Survival Stage (Stage 7)
    if args.stage in ["survival", "all"]:
        with log_stage(logger, "Survival Analysis & Prognostic Modeling"):
            run_survival_pipeline(args.config)

if __name__ == "__main__":
    main()
