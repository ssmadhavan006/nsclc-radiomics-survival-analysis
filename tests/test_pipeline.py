import os
import unittest
import numpy as np
import pandas as pd
import joblib
from src.utils import load_config
from src.model_utils import get_friendly_feature_name, format_ordinal, format_risk_markdown

class TestPipelineAndModel(unittest.TestCase):
    
    def test_config_loading(self):
        """Test that config file loads correctly and contains required fields."""
        config_path = "src/config.yaml"
        self.assertTrue(os.path.exists(config_path), "config.yaml not found")
        config = load_config(config_path)
        
        # Verify required configuration sections
        self.assertIn("paths", config)
        self.assertIn("preprocessing", config)
        self.assertIn("radiomics", config)
        self.assertIn("analysis", config)
        self.assertIn("survival", config)
        
        # Verify specific key paths
        self.assertIn("dataset_root", config["paths"])
        self.assertIn("clinical_csv", config["paths"])
        self.assertIn("output_dir", config["paths"])

    def test_friendly_name_mapping(self):
        """Test feature friendly name parsing and dictionary lookup."""
        # Test exact dictionary matches
        self.assertEqual(get_friendly_feature_name("original_shape_Elongation"), "Tumor Elongation")
        self.assertEqual(get_friendly_feature_name("wavelet-LLH_firstorder_Maximum"), "Maximum Density (LLH)")
        self.assertEqual(get_friendly_feature_name("wavelet-LLL_firstorder_Range"), "Density Range (LLL)")
        
        # Test fallback parser
        self.assertEqual(get_friendly_feature_name("wavelet-HLL_firstorder_Energy"), "Energy (firstorder HLL)")
        self.assertEqual(get_friendly_feature_name("custom_feature"), "feature (custom)")
        self.assertEqual(get_friendly_feature_name("custom"), "custom")

    def test_ordinal_formatting(self):
        """Test ordinal formatter output."""
        self.assertEqual(format_ordinal(1), "1st percentile")
        self.assertEqual(format_ordinal(22), "22nd percentile")
        self.assertEqual(format_ordinal(83), "83rd percentile")
        self.assertEqual(format_ordinal(11), "11th percentile")
        self.assertEqual(format_ordinal(12), "12th percentile")
        self.assertEqual(format_ordinal(13), "13th percentile")

    def test_risk_markdown_builder(self):
        """Test HTML styling of risk markdown builder."""
        risk_html = format_risk_markdown("High Risk", "+25% above median", "85th percentile")
        self.assertIn("High Risk", risk_html)
        self.assertIn("#f87171", risk_html) # Red color code for high risk
        self.assertIn("+25% above median", risk_html)
        
        risk_html_low = format_risk_markdown("Low Risk", "-10% below median", "34th percentile")
        self.assertIn("Low Risk", risk_html_low)
        self.assertIn("#4ade80", risk_html_low) # Green color code for low risk

    def test_serialized_checkpoints(self):
        """Verify existence and capability to load model checkpoints if they are generated."""
        scaler_path = "outputs/features/scaler.joblib"
        model_rad_path = "outputs/features/model_radiomics.joblib"
        model_comb_path = "outputs/features/model_combined.joblib"
        
        # This test passes if files are not generated yet (e.g. clean workspace) or loads and validates them
        if os.path.exists(scaler_path) and os.path.exists(model_rad_path) and os.path.exists(model_comb_path):
            scaler = joblib.load(scaler_path)
            model_rad = joblib.load(model_rad_path)
            model_comb = joblib.load(model_comb_path)
            
            # Check basic attributes
            self.assertTrue(hasattr(scaler, "transform"), "Scaler is invalid")
            self.assertTrue(hasattr(model_rad, "predict_partial_hazard"), "Model radiomics is invalid")
            self.assertTrue(hasattr(model_comb, "predict_partial_hazard"), "Model combined is invalid")

if __name__ == "__main__":
    unittest.main()
