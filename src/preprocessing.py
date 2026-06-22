import os
import glob
import pydicom
import SimpleITK as sitk
import numpy as np
import logging
from typing import Tuple, List, Optional, Dict, Any

logger = logging.getLogger("radiomics_pipeline")

def load_sorted_ct_series(ct_dir: str) -> Tuple[sitk.Image, List[float]]:
    """
    Reads a CT DICOM series, explicitly sorting files by ImagePositionPatient Z coordinate.
    
    Args:
        ct_dir: Path to the CT DICOM folder.
        
    Returns:
        Tuple: (sitk_ct_image, list_of_sorted_z_coordinates)
    """
    dcm_files = glob.glob(os.path.join(ct_dir, "*.dcm"))
    if not dcm_files:
        raise FileNotFoundError(f"No DICOM files found in CT directory: {ct_dir}")
        
    # Read each file header to get the Z coordinate for sorting (Rule 8)
    slice_positions = []
    for f in dcm_files:
        try:
            header = pydicom.dcmread(f, stop_before_pixels=True)
            ipp = header.ImagePositionPatient
            z_coord = float(ipp[2])
            slice_positions.append((f, z_coord))
        except Exception as e:
            raise ValueError(f"Failed to parse coordinate metadata from DICOM {f}: {str(e)}")
            
    # Sort files ascending by Z coordinate
    slice_positions.sort(key=lambda x: x[1])
    sorted_files = [x[0] for x in slice_positions]
    sorted_z = [x[1] for x in slice_positions]
    
    # Load using SimpleITK
    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(sorted_files)
    ct_image = reader.Execute()
    
    return ct_image, sorted_z

def load_and_align_dicom_seg(
    ct_image: sitk.Image, 
    sorted_z: List[float], 
    seg_path: str,
    target_segment_name: str = "Neoplasm, Primary"
) -> sitk.Image:
    """
    Loads a multi-frame DICOM SEG file and extracts/aligns the target segment 
    with the CT image using Coordinate Matching (Rule 9).
    
    Args:
        ct_image: The loaded and sorted SimpleITK CT image.
        sorted_z: Sorted list of Z-coordinates for the CT slices.
        seg_path: Path to the DICOM SEG file.
        target_segment_name: Name of the segment to extract.
        
    Returns:
        A SimpleITK binary mask image aligned exactly with ct_image.
    """
    dcm = pydicom.dcmread(seg_path)
    
    # Find the target segment index (1-based ReferencedSegmentNumber)
    target_seg_num = 1
    found_seg = False
    if hasattr(dcm, 'SegmentSequence'):
        for idx, seg in enumerate(dcm.SegmentSequence):
            label = getattr(seg, 'SegmentLabel', '')
            if target_segment_name.lower() in label.lower() or "primary" in label.lower() or "neoplasm" in label.lower():
                target_seg_num = idx + 1
                found_seg = True
                logger.info(f"Target segment '{target_segment_name}' found at index {target_seg_num} with label '{label}'")
                break
        if not found_seg:
            logger.warning(f"Target segment '{target_segment_name}' not found. Defaulting to Segment 1.")
            
    pixel_array = dcm.pixel_array
    frames = dcm.PerFrameFunctionalGroupsSequence
    
    # Map Z coordinate to the binary 2D frame mask
    z_to_frame_mask = {}
    for f_idx, frame in enumerate(frames):
        ref_seg = frame.SegmentIdentificationSequence[0].ReferencedSegmentNumber
        if ref_seg != target_seg_num:
            continue
            
        pps = frame.PlanePositionSequence[0]
        ipp = pps.ImagePositionPatient
        z_coord = float(ipp[2])
        z_to_frame_mask[round(z_coord, 3)] = pixel_array[f_idx]
        
    # Reconstruct the 3D binary mask array matched to CT slices Z positions
    mask_slices = []
    for z in sorted_z:
        z_rounded = round(z, 3)
        if z_rounded in z_to_frame_mask:
            mask_slices.append(z_to_frame_mask[z_rounded])
        else:
            # Empty slice if no segmentation exists for this coordinate
            mask_slices.append(np.zeros(ct_image.GetSize()[:2][::-1], dtype=np.uint8))
            
    mask_arr = np.stack(mask_slices, axis=0)
    # Ensure binary mask values (0 or 1)
    mask_arr = (mask_arr > 0).astype(np.uint8)
    
    # Convert numpy array to SimpleITK Image
    # NumPy order: (Z, Y, X) -> SimpleITK order: (X, Y, Z)
    mask_image = sitk.GetImageFromArray(mask_arr)
    mask_image.CopyInformation(ct_image)
    
    return mask_image

def resample_image(
    image: sitk.Image, 
    new_spacing: List[float], 
    interpolator: int
) -> sitk.Image:
    """
    Resamples a SimpleITK image to a target voxel spacing.
    
    Args:
        image: SimpleITK image to resample.
        new_spacing: Desired voxel spacing [dx, dy, dz] in mm.
        interpolator: SimpleITK interpolator choice (e.g. sitkBSpline, sitkNearestNeighbor).
        
    Returns:
        Resampled SimpleITK image.
    """
    original_spacing = image.GetSpacing()
    original_size = image.GetSize()
    
    if original_spacing == tuple(new_spacing):
        return image
        
    new_size = [
        int(round(original_size[i] * original_spacing[i] / new_spacing[i]))
        for i in range(3)
    ]
    
    resample = sitk.ResampleImageFilter()
    resample.SetInterpolator(interpolator)
    resample.SetOutputSpacing(new_spacing)
    resample.SetSize(new_size)
    resample.SetOutputDirection(image.GetDirection())
    resample.SetOutputOrigin(image.GetOrigin())
    resample.SetTransform(sitk.Transform())
    
    # Set default pixel value (especially important for HU values)
    if interpolator == sitk.sitkNearestNeighbor:
        resample.SetDefaultPixelValue(0)
    else:
        resample.SetDefaultPixelValue(-1000) # air HU
        
    return resample.Execute(image)

def preprocess_case(
    ct_dir: str, 
    seg_path: str, 
    config: Dict[str, Any]
) -> Tuple[sitk.Image, sitk.Image]:
    """
    Performs data loading, coordinate alignment, intensity clipping, and resampling (Stage 3).
    Verifies that spatial properties match exactly.
    
    Args:
        ct_dir: Path to CT DICOM directory.
        seg_path: Path to DICOM SEG file.
        config: Pipeline configuration dictionary.
        
    Returns:
        Tuple: (preprocessed_ct, preprocessed_mask)
    """
    # 1. Load CT series (sorted by Z position)
    ct_image, sorted_z = load_sorted_ct_series(ct_dir)
    
    # 2. Load and align Segmentation mask
    mask_image = load_and_align_dicom_seg(ct_image, sorted_z, seg_path)
    
    # 3. Clip CT intensity (HU Clipping) (Rule 12)
    hu_min = config["preprocessing"]["hu_min"]
    hu_max = config["preprocessing"]["hu_max"]
    clipped_ct = sitk.Clamp(ct_image, sitk.sitkFloat32, hu_min, hu_max)
    
    # 4. Resample CT and Mask to isotropic grid (Rule 12)
    target_spacing = config["preprocessing"]["target_spacing"]
    
    ct_interpolator_str = config["preprocessing"]["interpolator"]
    mask_interpolator_str = config["preprocessing"]["mask_interpolator"]
    
    # Convert string to SimpleITK interpolator
    ct_interpolator = getattr(sitk, ct_interpolator_str, sitk.sitkBSpline)
    mask_interpolator = getattr(sitk, mask_interpolator_str, sitk.sitkNearestNeighbor)
    
    resampled_ct = resample_image(clipped_ct, target_spacing, ct_interpolator)
    resampled_mask = resample_image(mask_image, target_spacing, mask_interpolator)
    
    # Ensure mask remains strictly binary (0 or 1)
    resampled_mask = sitk.BinaryThreshold(resampled_mask, 1, 255, 1, 0)
    
    # 5. Alignment check (Rule 9)
    spacing_tolerance = 1e-5
    origin_tolerance = 1e-3
    
    if resampled_ct.GetSize() != resampled_mask.GetSize():
        raise ValueError(f"Resampled CT and Mask size mismatch: CT {resampled_ct.GetSize()} vs Mask {resampled_mask.GetSize()}")
        
    if not all(abs(a - b) < spacing_tolerance for a, b in zip(resampled_ct.GetSpacing(), resampled_mask.GetSpacing())):
        raise ValueError(f"Resampled CT and Mask spacing mismatch: CT {resampled_ct.GetSpacing()} vs Mask {resampled_mask.GetSpacing()}")
        
    if not all(abs(a - b) < origin_tolerance for a, b in zip(resampled_ct.GetOrigin(), resampled_mask.GetOrigin())):
        raise ValueError(f"Resampled CT and Mask origin mismatch: CT {resampled_ct.GetOrigin()} vs Mask {resampled_mask.GetOrigin()}")
        
    if not all(abs(a - b) < 1e-4 for a, b in zip(resampled_ct.GetDirection(), resampled_mask.GetDirection())):
        raise ValueError("Resampled CT and Mask direction cosines mismatch.")
        
    # 6. Tumor mask overlap validation (Rule 10)
    # Check that mask has some active voxels (tumor region exists)
    statistics = sitk.LabelShapeStatisticsImageFilter()
    statistics.Execute(resampled_mask)
    if not statistics.HasLabel(1):
        raise ValueError("Resampled mask contains zero tumor voxels (label 1 is missing).")
        
    tumor_voxels = statistics.GetNumberOfPixels(1)
    min_voxels = config["preprocessing"]["min_tumor_voxels"]
    if tumor_voxels < min_voxels:
        raise ValueError(f"Tumor mask contains only {tumor_voxels} voxels, which is below the minimum threshold of {min_voxels}.")
        
    tumor_volume = statistics.GetPhysicalSize(1) # volume in mm³
    logger.info(f"Tumor mask validated: {tumor_voxels} voxels, volume: {tumor_volume:.2f} mm³")
    
    return resampled_ct, resampled_mask
