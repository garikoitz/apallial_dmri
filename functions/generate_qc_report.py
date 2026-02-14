#!/usr/bin/env python3
"""
Automated QC Report Generator for Tractography.

Iterates through all subjects in config1.yaml, finds generated tractography files,
creates screenshots using mrview (if not present), and compiles them into an HTML report.

Usage:
    python functions/generate_qc_report.py --group APTR
    python functions/generate_qc_report.py --all --force
    python functions/generate_qc_report.py --subject sub-Apallial17NmF

Author: MBRAIN Project
Pipeline: flip1
"""

import os
import argparse
import sys
import glob
import subprocess
import logging
from pathlib import Path
from datetime import datetime

# Set pipeline to flip1
os.environ['MBRAIN_PIPELINE'] = 'flip1'

# Add parent dir to path to import utils
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logger import setup_logger
from functions.config_loader import load_config, get_subject_groups

# Global logger
logger = logging.getLogger(__name__)


def get_tck_files(subject_dir):
    """Find all tract files in subject's segment_tracts directory."""
    tract_dir = os.path.join(subject_dir, "dwi", "segment_tracts")
    if not os.path.exists(tract_dir):
        return []
    
    # Priority: Cleaned > Raw
    # But we might want to see both? Let's just list all unique base names
    files = sorted(glob.glob(os.path.join(tract_dir, "*.tck")))
    return files


def generate_screenshot(tck_path, bg_path, output_png, force=False):
    """Generate screenshot using mrview."""
    if os.path.exists(output_png) and not force:
        return True
    
    # Create directory if needed
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    
    # mrview command
    # -mode 2: orthographic view
    # -capture.grab: take screenshot
    # -exit: close window immediately
    cmd = [
        "mrview", bg_path,
        "-tractography.load", tck_path,
        "-mode", "2",  # Orthographic
        "-capture.folder", os.path.dirname(output_png),
        "-capture.prefix", os.path.basename(output_png).replace(".png", ""),
        "-capture.grab",
        "-exit"
    ]
    
    # Note: On headless systems this might fail without xvfb-run
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # mrview capture filenames might have a number appended (e.g. prefix0000.png)
        # We need to rename it to exactly output_png if mrview added suffix
        base_prefix = os.path.basename(output_png).replace(".png", "")
        generated = glob.glob(os.path.join(os.path.dirname(output_png), f"{base_prefix}*.png"))
        
        if generated:
            # Sort to get the most recent/likely one
            generated.sort()
            actual_file = generated[-1]
            if actual_file != output_png:
                os.rename(actual_file, output_png)
            return True
        return False
        
    except subprocess.CalledProcessError:
        logger.error(f"Failed to generate screenshot for {os.path.basename(tck_path)}")
        return False


def generate_html_report(subjects_data, output_html):
    """Generate HTML report."""
    
    html = [
        "<!DOCTYPE html>",
        "<html><head>",
        "<title>Tractography QC Report</title>",
        "<style>",
        "body { font-family: sans-serif; margin: 20px; background: #f5f5f5; }",
        ".group-header { background: #333; color: white; padding: 10px; margin-top: 20px; }",
        ".subject-header { background: #666; color: white; padding: 5px 10px; margin-top: 10px; }",
        ".gallery { display: flex; flex-wrap: wrap; gap: 10px; }",
        ".card { background: white; padding: 10px; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); width: 300px; }",
        ".card img { width: 100%; height: auto; cursor: pointer; }",
        ".card h4 { margin: 5px 0; font-size: 12px; word-break: break-all; }",
        ".meta { font-size: 10px; color: #888; }",
        "</style>",
        "</head><body>",
        f"<h1>Tractography QC Report - {datetime.now().strftime('%Y-%m-%d %H:%M')}</h1>"
    ]
    
    for group, subjects in subjects_data.items():
        if not subjects: continue
        html.append(f"<h2 class='group-header'>{group}</h2>")
        
        for sub_id, tracts in subjects.items():
            if not tracts: continue
            html.append(f"<div class='subject-block'>")
            html.append(f"<h3 class='subject-header'>{sub_id} ({len(tracts)} tracts)</h3>")
            html.append("<div class='gallery'>")
            
            for tck_file, img_path in tracts:
                rel_img = os.path.relpath(img_path, os.path.dirname(output_html))
                name = os.path.basename(tck_file)
                html.append(f"""
                <div class='card'>
                    <a href='{rel_img}' target='_blank'>
                        <img src='{rel_img}' alt='{name}'>
                    </a>
                    <h4>{name}</h4>
                </div>
                """)
            
            html.append("</div></div>")
    
    html.append("</body></html>")
    
    with open(output_html, "w") as f:
        f.write("\n".join(html))
    
    return output_html


def main():
    parser = argparse.ArgumentParser(description="Generate QC report for tractography")
    parser.add_argument("--group", choices=["APTR", "AP", "CTRL"], help="Process specific group")
    parser.add_argument("--all", action="store_true", help="Process all groups")
    parser.add_argument("--subject", help="Process specific subject")
    parser.add_argument("--force", action="store_true", help="Force regenerate screenshots")
    parser.add_argument("--log-file", help="Path to log file")
    
    args = parser.parse_args()
    
    # Setup Logger
    global logger
    logger = setup_logger(__name__, args.log_file)
    
    if not (args.group or args.all or args.subject):
        parser.error("Must specify --group, --all, or --subject")
    
    config = load_config()
    groups_config = get_subject_groups()
    base_dir = config.get('paths', {}).get('derivatives_dir')
    repo_dir = config.get('paths', {}).get('repo_dir', '.')
    
    if not base_dir:
        logger.error("Could not find derivatives_dir in config")
        return

    # Determine subjects to process
    subjects_to_process = {} # group -> {sub_id: [(tck, img)]}
    
    target_groups = ["APTR", "AP", "CTRL"] if args.all else ([args.group] if args.group else [])
    
    # If specific subject, find its group
    if args.subject:
        found = False
        for g, subs in groups_config.items():
            for s in subs:
                if s['id'] == args.subject:
                    target_groups = [g]
                    subjects_to_process[g] = {s['id']: s['session']}
                    found = True
                    break
            if found: break
        if not found:
            logger.error(f"Subject {args.subject} not found")
            return
    else:
        # Load all subjects for target groups
        for g in target_groups:
            subjects_to_process[g] = {s['id']: s['session'] for s in groups_config.get(g, [])}

    # Output directory for QC
    qc_dir = os.path.join(repo_dir, "reports", "qc", "tractography")
    os.makedirs(qc_dir, exist_ok=True)
    
    # Data for HTML report
    report_data = {g: {} for g in target_groups}
    
    logger.info("Starting QC generation...")
    
    for group in target_groups:
        subs = subjects_to_process.get(group, {})
        for sub_id, ses_id in subs.items():
            logger.info(f"Processing {sub_id}...")
            
            sub_dir = os.path.join(base_dir, sub_id, ses_id)
            tck_files = get_tck_files(sub_dir)
            
            if not tck_files:
                logger.warning(f"  No tract files found for {sub_id}")
                continue
            
            # Find background image
            # Try T2w first, then FOD
            bg_path = os.path.join(sub_dir, "dwi", "t2w_highres_flip.mif")
            if not os.path.exists(bg_path):
                bg_path = os.path.join(sub_dir, "dwi", "wm_fod_norm.mif")
            
            if not os.path.exists(bg_path):
                logger.warning(f"  No background image found for {sub_id}")
                continue
                
            report_data[group][sub_id] = []
            
            # Sub-directory for images to keep things clean
            img_dir = os.path.join(qc_dir, "images", sub_id)
            os.makedirs(img_dir, exist_ok=True)
            
            for tck in tck_files:
                name = os.path.basename(tck).replace('.tck', '')
                img_path = os.path.join(img_dir, f"{name}.png")
                
                # Filter suggestions: Maybe skip SIFT2 weights files if they are just text?
                # No, get_tck_files only gets .tck
                
                if generate_screenshot(tck, bg_path, img_path, force=args.force):
                    report_data[group][sub_id].append((tck, img_path))
                    logger.info(f"  Generated: {name}")
                else:
                    if os.path.exists(img_path):
                        report_data[group][sub_id].append((tck, img_path))
                        logger.info(f"  Exists: {name}")
    
    # Generate HTML
    html_path = os.path.join(qc_dir, "index.html")
    generate_html_report(report_data, html_path)
    logger.info(f"QC Report generated at: {html_path}")

if __name__ == "__main__":
    main()
