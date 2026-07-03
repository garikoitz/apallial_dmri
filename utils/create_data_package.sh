#!/bin/bash
# create_data_package.sh — Build three OSF-ready ZIP components
#
# Components:
#   local -> atlas + manual ROI files
#   pre   -> BIDS + pre-mask mrtrix_newmask files, including edited masks
#   post  -> remaining mrtrix_newmask outputs
#
# Usage:
#   bash utils/create_data_package.sh
#   bash utils/create_data_package.sh --components local pre post
#   bash utils/create_data_package.sh --components local
#   bash utils/create_data_package.sh --components pre
#   bash utils/create_data_package.sh --components post
#   bash utils/create_data_package.sh --verify-only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
DATA_BASE="/export/home/glerma/glerma/PROYECTOS/mbrain"
STAGE_ROOT="/tmp/mbrain_package_staging"
OUTPUT_DIR="/tmp/mbrain_component_zips"
SELECTED_COMPONENTS=(local pre post)
VERIFY_ONLY=false

log_info() { echo -e "\033[0;34mℹ️  $1\033[0m"; }
log_ok() { echo -e "\033[0;32m✅ $1\033[0m"; }
log_warn() { echo -e "\033[1;33m⚠️  $1\033[0m"; }
log_error() { echo -e "\033[0;31m❌ $1\033[0m"; }

component_enabled() {
    local wanted="$1"
    for component in "${SELECTED_COMPONENTS[@]}"; do
        [[ "$component" == "$wanted" ]] && return 0
    done
    return 1
}

zip_stage() {
    local component="$1"
    local stage_dir="$2"
    local zip_path="$3"

    rm -rf "$stage_dir"
    mkdir -p "$stage_dir"

    case "$component" in
        local)
            mkdir -p "$stage_dir/local"
            cp -a "$DATA_BASE/local/atlas" "$stage_dir/local/"
            cp -a "$DATA_BASE/local/manual_roi_Slicer" "$stage_dir/local/"
            cp -a "$DATA_BASE/local/manual_roi_flip" "$stage_dir/local/"
            printf 'local\n' > "$stage_dir/.mbrain_component_local"
            ;;
        pre)
            mkdir -p "$stage_dir/data/bids"
            cp -a "$DATA_BASE/data/bids/." "$stage_dir/data/bids/"

            mkdir -p "$stage_dir/data/derivatives/mrtrix_newmask"
            rsync -a --prune-empty-dirs \
                --include='*/' \
                --include='t2w_highres_flip*' \
                --include='dwi*' \
                --include='mask*' \
                --include='*.bval' \
                --include='*.bvec' \
                --include='*.json' \
                --exclude='*' \
                "$DATA_BASE/data/derivatives/mrtrix_newmask/" \
                "$stage_dir/data/derivatives/mrtrix_newmask/"

            printf 'pre\n' > "$stage_dir/.mbrain_component_pre"
            ;;
        post)
            mkdir -p "$stage_dir/data/derivatives/mrtrix_newmask"
            rsync -a --prune-empty-dirs \
                --include='*/' \
                --exclude='t2w_highres_flip*' \
                --exclude='dwi*' \
                --exclude='mask*' \
                --exclude='*.bval' \
                --exclude='*.bvec' \
                "$DATA_BASE/data/derivatives/mrtrix_newmask/" \
                "$stage_dir/data/derivatives/mrtrix_newmask/"

            printf 'post\n' > "$stage_dir/.mbrain_component_post"
            ;;
        *)
            log_error "Unknown component: $component"
            exit 1
            ;;
    esac

    ( cd "$stage_dir" && zip -r -9 -q "$zip_path" . )

    if unzip -t "$zip_path" >/dev/null 2>&1; then
        log_ok "ZIP integrity verified for $component"
    else
        log_error "ZIP integrity check failed for $component"
        exit 1
    fi

    sha256sum "$zip_path" > "${zip_path}.sha256"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --verify-only)
            VERIFY_ONLY=true
            shift
            ;;
        --components)
            shift
            SELECTED_COMPONENTS=()
            while [[ $# -gt 0 ]]; do
                case "$1" in
                    local|pre|post)
                        SELECTED_COMPONENTS+=("$1")
                        shift
                        ;;
                    all)
                        SELECTED_COMPONENTS=(local pre post)
                        shift
                        ;;
                    --*)
                        break
                        ;;
                    *)
                        break
                        ;;
                esac
            done
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║  MBRAIN Component Packaging for OSF.io                            ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo

log_info "Repository: $REPO_DIR"
log_info "Data base:  $DATA_BASE"
log_info "Components: ${SELECTED_COMPONENTS[*]}"

if [[ ! -d "$DATA_BASE" ]]; then
    log_error "Data base not found: $DATA_BASE"
    exit 1
fi

required_checks=(
    "$DATA_BASE/local/atlas/average_template_25.nii.gz:Allen CCF template"
    "$DATA_BASE/local/atlas/annotation_25.nii.gz:Allen CCF annotation"
    "$DATA_BASE/local/manual_roi_Slicer:Manual ROI directory"
    "$DATA_BASE/local/manual_roi_flip:Flipped ROI directory"
    "$DATA_BASE/data/bids:Raw BIDS directory"
    "$DATA_BASE/data/derivatives/mrtrix_newmask:Newmask derivatives"
)

for spec in "${required_checks[@]}"; do
    path="${spec%%:*}"
    label="${spec##*:}"
    if [[ -e "$path" ]]; then
        log_ok "$label"
    else
        log_error "$label not found: $path"
        exit 1
    fi
done

if [[ "$VERIFY_ONLY" == true ]]; then
    log_ok "Verification passed"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"

if component_enabled local; then
    log_info "Packaging local component"
    zip_stage local "$STAGE_ROOT/local" "$OUTPUT_DIR/mbrain_local.zip"
fi

if component_enabled pre; then
    log_info "Packaging pre-mask component"
    zip_stage pre "$STAGE_ROOT/pre" "$OUTPUT_DIR/mbrain_pre.zip"
fi

if component_enabled post; then
    log_info "Packaging post-mask component"
    zip_stage post "$STAGE_ROOT/post" "$OUTPUT_DIR/mbrain_post.zip"
fi

echo
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║  PACKAGING COMPLETE                                                ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo

for zip_path in "$OUTPUT_DIR"/*.zip; do
    [[ -e "$zip_path" ]] || continue
    size=$(du -sh "$zip_path" | awk '{print $1}')
    echo "  $(basename "$zip_path"): $size"
    echo "  $(basename "$zip_path").sha256"
done
echo
echo "Upload the three ZIPs to OSF and install them with utils/setup_data.py --components"
DATA_BASE="/export/home/glerma/glerma/PROYECTOS/mbrain"
STAGE_ROOT="/tmp/mbrain_package_staging"
OUTPUT_DIR="/tmp/mbrain_component_zips"
VERIFY_ONLY=false
SELECTED_COMPONENTS=(local pre post)

log_info() { echo -e "\033[0;34mℹ️  $1\033[0m"; }
log_ok() { echo -e "\033[0;32m✅ $1\033[0m"; }
log_warn() { echo -e "\033[1;33m⚠️  $1\033[0m"; }
log_error() { echo -e "\033[0;31m❌ $1\033[0m"; }

component_enabled() {
    local wanted="$1"
    for comp in "${SELECTED_COMPONENTS[@]}"; do
        [[ "$comp" == "$wanted" ]] && return 0
    done
    return 1
}

copy_local_component() {
    local stage_dir="$1"
    mkdir -p "$stage_dir/local"
    cp -a "$DATA_BASE/local/atlas" "$stage_dir/local/"
    cp -a "$DATA_BASE/local/manual_roi_Slicer" "$stage_dir/local/"
    cp -a "$DATA_BASE/local/manual_roi_flip" "$stage_dir/local/"
}

copy_pre_component() {
    local stage_dir="$1"
    mkdir -p "$stage_dir/data/bids"
    cp -a "$DATA_BASE/data/bids/." "$stage_dir/data/bids/"

    local source_root="$DATA_BASE/data/derivatives/mrtrix_newmask"
    local target_root="$stage_dir/data/derivatives/mrtrix_newmask"
    mkdir -p "$target_root"

    find "$source_root" -type f \( \
        -name 't2w_highres_flip.mif' -o \
        -name 't2w_highres_flip.nii.gz' -o \
        -name 'dwi.mif' -o \
        -name 'dwi.nii.gz' -o \
        -name 'dwi*.mif' -o \
        -name 'dwi*.nii.gz' -o \
        -name '*.bval' -o \
        -name '*.bvec' -o \
        -name 'mask.mif' -o \
        -name 'mask.nii.gz' -o \
        -name 'mask_anat.mif' -o \
        -name 'mask_anat.nii.gz' -o \
        -name 'mask_new.mif' -o \
        -name 'mask_new.nii.gz' -o \
        -name 'mask_new_edited.mif' -o \
        -name 'mask_new_edited.nii.gz' \
    \) | while IFS= read -r file; do
        rel="${file#${source_root}/}"
        mkdir -p "$(dirname "$target_root/$rel")"
        cp "$file" "$target_root/$rel"
    done
}

copy_post_component() {
    local stage_dir="$1"
    local source_root="$DATA_BASE/data/derivatives/mrtrix_newmask"
    local target_root="$stage_dir/data/derivatives/mrtrix_newmask"
    mkdir -p "$target_root"

    find "$source_root" -type f \
        ! -name 't2w_highres_flip.mif' \
        ! -name 't2w_highres_flip.nii.gz' \
        ! -name 'dwi.mif' \
        ! -name 'dwi.nii.gz' \
        ! -name 'dwi*.mif' \
        ! -name 'dwi*.nii.gz' \
        ! -name '*.bval' \
        ! -name '*.bvec' \
        ! -name 'mask.mif' \
        ! -name 'mask.nii.gz' \
        ! -name 'mask_anat.mif' \
        ! -name 'mask_anat.nii.gz' \
        ! -name 'mask_new.mif' \
        ! -name 'mask_new.nii.gz' \
        ! -name 'mask_new_edited.mif' \
        ! -name 'mask_new_edited.nii.gz' | while IFS= read -r file; do
        rel="${file#${source_root}/}"
        mkdir -p "$(dirname "$target_root/$rel")"
        cp "$file" "$target_root/$rel"
    done
}

zip_component() {
    local component="$1"
    local stage_dir="$2"
    local zip_path="$3"

    rm -rf "$stage_dir"
    mkdir -p "$stage_dir"

    case "$component" in
        local) copy_local_component "$stage_dir" ;;
        pre) copy_pre_component "$stage_dir" ;;
        post) copy_post_component "$stage_dir" ;;
        *) log_error "Unknown component: $component"; exit 1 ;;
    esac

    ( cd "$stage_dir" && zip -r -9 -q "$zip_path" . )

    if unzip -t "$zip_path" >/dev/null 2>&1; then
        log_ok "ZIP integrity verified for $component"
    else
        log_error "ZIP integrity check failed for $component"
        exit 1
    fi

    sha256sum "$zip_path" > "${zip_path}.sha256"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --verify-only)
            VERIFY_ONLY=true
            shift
            ;;
        --components)
            shift
            SELECTED_COMPONENTS=()
            while [[ $# -gt 0 ]]; do
                case "$1" in
                    local|pre|post)
                        SELECTED_COMPONENTS+=("$1")
                        shift
                        ;;
                    all)
                        SELECTED_COMPONENTS=(local pre post)
                        shift
                        ;;
                    --*) break ;;
                    *) break ;;
                esac
            done
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║  MBRAIN Component Packaging for OSF.io                            ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo
log_info "Repository: $REPO_DIR"
log_info "Data base:  $DATA_BASE"
log_info "Components: ${SELECTED_COMPONENTS[*]}"

if [[ ! -d "$DATA_BASE" ]]; then
    log_error "Data base not found: $DATA_BASE"
    exit 1
fi

required_checks=(
    "$DATA_BASE/local/atlas/average_template_25.nii.gz:Allen CCF template"
    "$DATA_BASE/local/atlas/annotation_25.nii.gz:Allen CCF annotation"
    "$DATA_BASE/local/manual_roi_Slicer:Manual ROI directory"
    "$DATA_BASE/local/manual_roi_flip:Flipped ROI directory"
    "$DATA_BASE/data/bids:Raw BIDS directory"
    "$DATA_BASE/data/derivatives/mrtrix_newmask:Newmask derivatives"
)

for spec in "${required_checks[@]}"; do
    path="${spec%%:*}"
    label="${spec##*:}"
    [[ -e "$path" ]] && log_ok "$label" || { log_error "$label not found: $path"; exit 1; }
done

if [[ "$VERIFY_ONLY" == true ]]; then
    log_ok "Verification passed"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"

if component_enabled local; then
    log_info "Packaging local component"
    zip_component local "$STAGE_ROOT/local" "$OUTPUT_DIR/mbrain_local.zip"
fi

if component_enabled pre; then
    log_info "Packaging pre-mask component"
    zip_component pre "$STAGE_ROOT/pre" "$OUTPUT_DIR/mbrain_pre.zip"
fi

if component_enabled post; then
    log_info "Packaging post-mask component"
    zip_component post "$STAGE_ROOT/post" "$OUTPUT_DIR/mbrain_post.zip"
fi

echo
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║  PACKAGING COMPLETE                                                ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo
for zip_path in "$OUTPUT_DIR"/*.zip; do
    [[ -e "$zip_path" ]] || continue
    size=$(du -sh "$zip_path" | awk '{print $1}')
    echo "  $(basename "$zip_path"): $size"
    echo "  $(basename "$zip_path").sha256"
done
echo
echo "Upload the three ZIPs to OSF and install selected components with utils/setup_data.py --components"
