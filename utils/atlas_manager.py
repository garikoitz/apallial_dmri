#!/usr/bin/env python3
"""
Atlas Manager for Mouse Brain Analysis

Downloads and manages the Allen Mouse Brain Common Coordinate Framework (CCFv3)
and provides utilities for atlas-based coregistration.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import subprocess
import requests

logger = logging.getLogger(__name__)


class AllenAtlasManager:
    """Manager for Allen Mouse Brain CCF atlas."""
    
    # Structure IDs for key regions (from Allen ontology)
    STRUCTURE_IDS = {
        'cortex': 315,  # Isocortex
        'striatum': 477,  # Striatum
        'thalamus': 549,  # Thalamus
        'hippocampus': 1089,  # Hippocampal formation
        'cerebellum': 512,  # Cerebellum
    }
    
    # Hemisphere suffixes for bilateral structures
    HEMISPHERES = ['left', 'right']
    
    def __init__(self, atlas_dir: Optional[Path] = None):
        """
        Initialize atlas manager.
        
        Parameters
        ----------
        atlas_dir : Path, optional
            Directory to store atlas files. If None, uses data/atlas/
        """
        if atlas_dir is None:
            from .config import Config
            atlas_dir = Config.DATA_DIR / 'atlas' / 'allen_ccf'
        
        self.atlas_dir = Path(atlas_dir)
        self.atlas_dir.mkdir(parents=True, exist_ok=True)
        
        # Atlas files at 25um resolution (good balance for mouse DTI)
        self.resolution = 25
        self.template_file = self.atlas_dir / f'average_template_{self.resolution}.nii.gz'
        self.annotation_file = self.atlas_dir / f'annotation_{self.resolution}.nii.gz'
        self.structures_file = self.atlas_dir / 'structures.json'
        
        # Base URL for Allen CCF downloads
        self.base_url = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf'
    
    def is_downloaded(self) -> bool:
        """Check if atlas is already downloaded."""
        return (self.template_file.exists() and 
                self.annotation_file.exists() and
                self.structures_file.exists())
    
    def download(self, force: bool = False) -> Dict[str, Path]:
        """
        Download Allen CCF atlas files.
        
        Parameters
        ----------
        force : bool
            Force re-download even if files exist
        
        Returns
        -------
        Dict[str, Path]
            Paths to downloaded files
        """
        if self.is_downloaded() and not force:
            logger.info(f"Atlas already downloaded at {self.atlas_dir}")
            return self.get_atlas_files()
        
        logger.info(f"Downloading Allen Mouse Brain CCF (resolution: {self.resolution}µm)")
        logger.info(f"Target directory: {self.atlas_dir}")
        
        # Download template
        logger.info("Downloading average template...")
        template_url = f"{self.base_url}/average_template/average_template_{self.resolution}.nrrd"
        self._download_and_convert(template_url, self.template_file)
        
        # Download annotation
        logger.info("Downloading annotation volume...")
        annotation_url = f"{self.base_url}/annotation/ccf_2017/annotation_{self.resolution}.nrrd"
        self._download_and_convert(annotation_url, self.annotation_file)
        
        # Download structure ontology
        logger.info("Downloading structure ontology...")
        structures_url = "http://api.brain-map.org/api/v2/structure_graph_download/1.json"
        self._download_json(structures_url, self.structures_file)
        
        logger.info("✓ Allen CCF download complete!")
        return self.get_atlas_files()
    
    def _download_and_convert(self, url: str, output_file: Path):
        """Download NRRD file and convert to NIfTI."""
        import nrrd
        import nibabel as nib
        
        nrrd_file = output_file.with_suffix('.nrrd')
        
        # Download
        logger.info(f"  Downloading from {url}")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        last_percent = 0
        with open(nrrd_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    percent = int((downloaded / total_size) * 100)
                    # Only log every 10%
                    if percent >= last_percent + 10:
                        logger.info(f"    Progress: {percent}%")
                        last_percent = percent
        
        logger.info(f"    Progress: 100%")
        
        # Convert to NIfTI using nrrd and nibabel
        logger.info(f"  Converting to NIfTI: {output_file.name}")
        data, header = nrrd.read(str(nrrd_file))
        
        # Create proper affine matrix from NRRD space directions
        # Allen CCF NRRD uses LPS orientation, NIfTI uses RAS+
        # Need to convert: LPS → RAS (flip X and Y)
        import numpy as np
        
        # Get space directions from NRRD header
        space_directions = header.get('space directions', None)
        space_origin = header.get('space origin', None)
        space = header.get('space', 'unknown')
        
        logger.info(f"  NRRD space: {space}")
        logger.info(f"  Original data shape: {data.shape}")
        
        # For Allen CCF: the data is stored as (left-right, posterior-anterior, superior-inferior)
        # But we need to check the actual axis order and reorient to RAS
        
        if space.lower() == 'left-posterior-superior':
            # The data needs to be reoriented from LPS to RAS
            # LPS array order: (L-R, P-A, S-I) where L=+, P=+, S=+
            # RAS needs: (R-L, A-P, S-I) where R=+, A=+, S=+
            
            # Flip the data along X and Y axes to convert LPS to RAS
            data = np.flip(data, axis=0)  # Flip L-R to R-L
            data = np.flip(data, axis=1)  # Flip P-A to A-P
            # Z axis (S-I) stays the same
            
            logger.info(f"  Flipped data to RAS orientation")
        
        # Create simple RAS+ affine with correct spacing
        if space_directions is not None:
            # Get voxel sizes (convert from µm to mm)
            voxel_sizes = np.array([
                np.linalg.norm(space_directions[0]),
                np.linalg.norm(space_directions[1]),
                np.linalg.norm(space_directions[2])
            ]) / 1000.0
        else:
            voxel_sizes = np.array([0.025, 0.025, 0.025])
        
        # Create affine matrix - simple diagonal for RAS orientation
        affine = np.diag([voxel_sizes[0], voxel_sizes[1], voxel_sizes[2], 1.0])
        
        # Center the image
        affine[0, 3] = -data.shape[0] * voxel_sizes[0] / 2
        affine[1, 3] = -data.shape[1] * voxel_sizes[1] / 2
        affine[2, 3] = -data.shape[2] * voxel_sizes[2] / 2
        
        logger.info(f"  Final data shape: {data.shape}")
        logger.info(f"  Voxel sizes: {voxel_sizes} mm")
        logger.info(f"  Affine matrix (RAS+):\n{affine}")
        
        # Create NIfTI image and save directly
        # NOTE: Subject data is flipped at preprocessing time to match atlas
        # Therefore NO flip is applied to the atlas itself
        nifti_img = nib.Nifti1Image(data, affine=affine)
        nib.save(nifti_img, str(output_file))
        
        # Clean up temporary NRRD file
        nrrd_file.unlink()
        
        logger.info(f"  ✓ Converted to {output_file.name} (LPS→RAS only, no additional flip)")

    
    def _download_json(self, url: str, output_file: Path):
        """Download JSON file."""
        logger.info(f"  Downloading from {url}")
        response = requests.get(url)
        response.raise_for_status()
        
        with open(output_file, 'w') as f:
            json.dump(response.json(), f, indent=2)
    
    def get_atlas_files(self) -> Dict[str, Path]:
        """Get paths to atlas files."""
        return {
            'template': self.template_file,
            'annotation': self.annotation_file,
            'structures': self.structures_file,
        }
    
    def load_structures(self) -> Dict:
        """Load structure ontology."""
        with open(self.structures_file, 'r') as f:
            return json.load(f)
    
    def get_structure_mask(self, structure_name: str, output_file: Path) -> Path:
        """
        Extract binary mask for a specific structure.
        
        Parameters
        ----------
        structure_name : str
            Name of structure ('cortex', 'striatum', 'thalamus', etc.)
        output_file : Path
            Output file path for mask
        
        Returns
        -------
        Path
            Path to created mask file
        """
        if structure_name not in self.STRUCTURE_IDS:
            raise ValueError(f"Unknown structure: {structure_name}. "
                           f"Available: {list(self.STRUCTURE_IDS.keys())}")
        
        structure_id = self.STRUCTURE_IDS[structure_name]
        
        logger.info(f"Extracting {structure_name} mask (ID: {structure_id})")
        
        # Load structure ontology to get all child IDs
        structures = self.load_structures()
        child_ids = self._get_all_child_ids(structures, structure_id)
        
        logger.info(f"  Found {len(child_ids)} sub-structures")
        
        # Create mask using numpy - much faster and more reliable than mrcalc for many IDs
        import nibabel as nib
        import numpy as np
        
        # Load annotation volume
        logger.info(f"  Loading annotation volume...")
        annotation_img = nib.load(self.annotation_file)
        annotation_data = annotation_img.get_fdata().astype(int)
        
        # Create binary mask
        logger.info(f"  Creating binary mask...")
        mask_data = np.isin(annotation_data, child_ids).astype(np.uint8)
        
        # Count non-zero voxels
        n_voxels = np.sum(mask_data > 0)
        logger.info(f"  Mask contains {n_voxels:,} voxels")
        
        if n_voxels == 0:
            logger.warning(f"  Warning: {structure_name} mask is empty! Check structure ID {structure_id}")
        
        # Save as NIfTI
        mask_img = nib.Nifti1Image(mask_data, annotation_img.affine, annotation_img.header)
        nib.save(mask_img, output_file)
        
        # Standardize strides to [1,2,3] using mrconvert
        temp_file = output_file.parent / f"temp_{output_file.name}"
        subprocess.run([
            'mrconvert',
            str(output_file),
            str(temp_file),
            '-strides', '1,2,3',
            '-force'
        ], check=True, capture_output=True)
        
        # Replace original with standardized version
        temp_file.replace(output_file)
        
        logger.info(f"  ✓ Mask saved to {output_file}")
        return output_file
    
    def split_hemisphere(self, mask_file: Path, output_left: Path, output_right: Path) -> Tuple[Path, Path]:
        """
        Split a bilateral mask into left and right hemispheres.
        
        The Allen CCF uses RAS+ orientation where the sagittal midline divides
        left (negative X) and right (positive X) hemispheres.
        
        Parameters
        ----------
        mask_file : Path
            Input bilateral mask
        output_left : Path
            Output file for left hemisphere mask
        output_right : Path
            Output file for right hemisphere mask
        
        Returns
        -------
        Tuple[Path, Path]
            Paths to left and right hemisphere masks
        """
        import nibabel as nib
        import numpy as np
        
        # Load mask
        img = nib.load(mask_file)
        data = img.get_fdata()
        
        # Get image center along X axis (sagittal)
        x_center = data.shape[0] // 2
        
        # Split into hemispheres
        left_data = np.zeros_like(data)
        right_data = np.zeros_like(data)
        
        # Left hemisphere: X < center
        left_data[:x_center, :, :] = data[:x_center, :, :]
        
        # Right hemisphere: X >= center
        right_data[x_center:, :, :] = data[x_center:, :, :]
        
        # Save hemispheres
        left_img = nib.Nifti1Image(left_data, img.affine, img.header)
        right_img = nib.Nifti1Image(right_data, img.affine, img.header)
        
        nib.save(left_img, output_left)
        nib.save(right_img, output_right)
        
        logger.info(f"  ✓ Left hemisphere saved to {output_left}")
        logger.info(f"  ✓ Right hemisphere saved to {output_right}")
        
        return output_left, output_right
    
    def _get_all_child_ids(self, structures: Dict, parent_id: int) -> List[int]:
        """
        Get all descendant structure IDs.
        
        Parameters
        ----------
        structures : Dict
            Structure ontology from Allen
        parent_id : int
            Parent structure ID
        
        Returns
        -------
        List[int]
            List of all child structure IDs (including parent)
        """
        child_ids = []
        
        def find_structure_and_descendants(struct, target_id):
            """Recursively find a structure by ID and get all its descendants."""
            if struct['id'] == target_id:
                # Found the target structure - add it and all descendants
                child_ids.append(struct['id'])
                if 'children' in struct and struct['children']:
                    for child in struct['children']:
                        add_all_descendants(child)
                return True
            
            # Keep searching in children
            if 'children' in struct and struct['children']:
                for child in struct['children']:
                    if find_structure_and_descendants(child, target_id):
                        return True
            return False
        
        def add_all_descendants(struct):
            """Add a structure and all its descendants to the list."""
            child_ids.append(struct['id'])
            if 'children' in struct and struct['children']:
                for child in struct['children']:
                    add_all_descendants(child)
        
        # Start from root of structure graph
        if 'msg' in structures:
            for root in structures['msg']:
                if find_structure_and_descendants(root, parent_id):
                    break
        
        return child_ids
    
    def extract_all_rois(self, output_dir: Path) -> Dict[str, Path]:
        """
        Extract all major ROIs (cortex, striatum, thalamus) split by hemisphere.
        
        ROIs are kept at full atlas resolution. ANTs will handle resolution 
        differences during registration automatically.
        
        Parameters
        ----------
        output_dir : Path
            Directory to save ROI masks
        
        Returns
        -------
        Dict[str, Path]
            Dictionary mapping structure names to mask file paths
            Keys: 'cortex_left', 'cortex_right', 'striatum_left', etc.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        rois = {}
        for structure_name in ['cortex', 'striatum', 'thalamus']:
            # Extract bilateral mask at full resolution
            bilateral_file = output_dir / f'atlas_{structure_name}_bilateral.nii.gz'
            self.get_structure_mask(structure_name, bilateral_file)
            
            # Split into hemispheres (keeping full resolution)
            left_file = output_dir / f'atlas_{structure_name}_left.nii.gz'
            right_file = output_dir / f'atlas_{structure_name}_right.nii.gz'
            
            logger.info(f"  Splitting {structure_name} into hemispheres...")
            self.split_hemisphere(bilateral_file, left_file, right_file)
            
            # Store both hemispheres
            rois[f'{structure_name}_left'] = left_file
            rois[f'{structure_name}_right'] = right_file
        
        return rois


def main():
    """Command-line interface for atlas management."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Download and manage Allen Mouse Brain CCF')
    parser.add_argument('--download', action='store_true',
                       help='Download atlas files')
    parser.add_argument('--extract-rois', action='store_true',
                       help='Extract cortex, striatum, and thalamus masks')
    parser.add_argument('--force', action='store_true',
                       help='Force re-download')
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s - %(levelname)s - %(message)s')
    
    manager = AllenAtlasManager()
    
    if args.download or not manager.is_downloaded():
        manager.download(force=args.force)
    
    if args.extract_rois:
        atlas_files = manager.get_atlas_files()
        roi_dir = manager.atlas_dir / 'rois'
        rois = manager.extract_all_rois(roi_dir)
        print(f"\nExtracted ROIs:")
        for name, path in rois.items():
            print(f"  {name}: {path}")


if __name__ == '__main__':
    main()
