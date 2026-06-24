import os
import time
import logging
import yaml
from typing import Any, Dict, Generator
from contextlib import contextmanager

def setup_logger(log_dir: str = "outputs/logs") -> logging.Logger:
    """
    Sets up the pipeline logger to write to both console and a timestamped file.
    
    Args:
        log_dir: The directory where logs will be stored.
        
    Returns:
        A configured logging.Logger instance.
    """
    os.makedirs(log_dir, exist_ok=True)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"pipeline_{timestamp}.log")
    
    logger = logging.getLogger("radiomics_pipeline")
    logger.setLevel(logging.INFO)
    
    # Avoid duplicate handlers if setup is called multiple times
    if not logger.handlers:
        # File handler
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_formatter = logging.Formatter(
            "[%(asctime)s] [%(levelname)s] [%(filename)s:%(lineno)d] - %(message)s"
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_formatter = logging.Formatter(
            "[%(asctime)s] [%(levelname)s] - %(message)s", datefmt="%H:%M:%S"
        )
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        
    return logger

def load_config(config_path: str = "src/config.yaml") -> Dict[str, Any]:
    """
    Loads and parses the master configuration YAML file.
    
    Args:
        config_path: Path to the config.yaml file.
        
    Returns:
        A dictionary containing configurations.
        
    Raises:
        FileNotFoundError: If the configuration file does not exist.
        yaml.YAMLError: If parsing configuration fails.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found at: {config_path}")
        
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        return config
    except yaml.YAMLError as e:
        logger = logging.getLogger("radiomics_pipeline")
        logger.error(f"Failed to parse config.yaml: {str(e)}")
        raise

@contextmanager
def log_stage(logger: logging.Logger, stage_name: str) -> Generator[None, None, None]:
    """
    A context manager to log the start, finish, duration, and status of a pipeline stage.
    
    Args:
        logger: Logger instance to use.
        stage_name: Name of the pipeline stage.
    """
    logger.info(f"=== Starting Stage: {stage_name} ===")
    start_time = time.time()
    try:
        yield
        duration = time.time() - start_time
        logger.info(f"=== Finished Stage: {stage_name} | Status: SUCCESS | Duration: {duration:.2f}s ===")
    except Exception as e:
        duration = time.time() - start_time
        logger.error(f"=== Failed Stage: {stage_name} | Status: ERROR | Duration: {duration:.2f}s ===")
        logger.error(f"Exception details: {str(e)}", exc_info=True)
        raise e
