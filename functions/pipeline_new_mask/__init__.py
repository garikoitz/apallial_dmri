"""Utilities for refactoring notebooks/pipeline_groups_new_mask.py into testable modules."""

from .exec_utils import check_file, launch_gui, run_command, run_parallel
from .mask_generation import (
    DEFAULT_GROUP_ATLAS_IDS,
    build_group_atlas_info,
    build_mask_generation_args,
    copy_atlas_mask,
    register_and_warp_mask,
)
from .fod_verify import collect_fod_verification, collect_group_response_status
from .mask_install import (
    build_mask_grid_report,
    build_mask_install_plan,
    compute_pipeline_checkpoint,
    resolve_dwi_template,
    select_new_mask,
)
from .rerun_tasks import (
    build_batch_process_tasks,
    build_ctrl_cleaning_tasks,
    build_ctrl_passthrough_extraction_tasks,
    build_cst_cleaning_tasks,
    build_cst_passthrough_extraction_tasks,
    build_cst_tckgen_tasks,
    build_passthrough_cleaning_tasks,
    build_passthrough_extraction_tasks,
    build_shrink_segments_command,
    build_sift2_tasks,
    build_tck_downsample_tasks,
    build_wholebrain_tractography_tasks,
    collect_passthrough_tract_counts,
    collect_ctrl_passthrough_counts,
)
from .roi_qc import (
    evaluate_laterality_consistency,
    get_grid_info,
    get_roi_pairs_for_group,
    load_roi_data,
    select_sample_tract,
    select_raw_cleaned_comparison_pair,
)
from .report_plots import draw_single_seg_panel
from .cst_roi import build_cst_regrid_plan, collect_cst_status, resolve_cst_source_rois
from .cst_roi import build_cst_final_qc_summary
from .cst_roi import collect_cst_tckgen_counts
from .cst_roi import build_temp_cst_exclusion_tasks, collect_temp_exclusion_rois
from .cst_roi import collect_cst_cleaned_counts, collect_cst_raw_tract_counts
from .trk_conversion import convert_tck_to_trk
from .subjects import build_subject_records, summarize_subjects

__all__ = [
    "run_command",
    "check_file",
    "launch_gui",
    "run_parallel",
    "DEFAULT_GROUP_ATLAS_IDS",
    "build_group_atlas_info",
    "build_mask_generation_args",
    "register_and_warp_mask",
    "copy_atlas_mask",
    "select_new_mask",
    "resolve_dwi_template",
    "build_mask_install_plan",
    "build_mask_grid_report",
    "compute_pipeline_checkpoint",
    "build_batch_process_tasks",
    "build_ctrl_passthrough_extraction_tasks",
    "build_ctrl_cleaning_tasks",
    "build_cst_cleaning_tasks",
    "build_cst_passthrough_extraction_tasks",
    "build_cst_tckgen_tasks",
    "build_passthrough_cleaning_tasks",
    "build_shrink_segments_command",
    "build_passthrough_extraction_tasks",
    "build_wholebrain_tractography_tasks",
    "build_sift2_tasks",
    "build_tck_downsample_tasks",
    "collect_passthrough_tract_counts",
    "collect_ctrl_passthrough_counts",
    "collect_fod_verification",
    "collect_group_response_status",
    "get_roi_pairs_for_group",
    "load_roi_data",
    "get_grid_info",
    "evaluate_laterality_consistency",
    "select_sample_tract",
    "select_raw_cleaned_comparison_pair",
    "draw_single_seg_panel",
    "resolve_cst_source_rois",
    "collect_cst_status",
    "build_cst_regrid_plan",
    "build_cst_final_qc_summary",
    "collect_cst_tckgen_counts",
    "collect_temp_exclusion_rois",
    "build_temp_cst_exclusion_tasks",
    "collect_cst_cleaned_counts",
    "collect_cst_raw_tract_counts",
    "convert_tck_to_trk",
    "build_subject_records",
    "summarize_subjects",
]
