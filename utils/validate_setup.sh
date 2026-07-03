#!/bin/bash
# validate_setup.sh — Quick validation that everything is ready
#
# Usage:
#   bash utils/validate_setup.sh

set -e

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}   MBRAIN Setup Validation${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}\n"

pass=0
fail=0

check() {
    local name="$1"
    local condition="$2"
    
    if eval "$condition"; then
        echo -e "${GREEN}✅${NC} $name"
        ((pass++))
    else
        echo -e "${RED}❌${NC} $name"
        ((fail++))
    fi
}

# Repository structure
echo -e "${BLUE}Repository Structure${NC}"
check "Repository root exists" "[[ -d '$REPO_DIR' ]]"
check "notebooks/ directory exists" "[[ -d '$REPO_DIR/notebooks' ]]"
check "functions/ directory exists" "[[ -d '$REPO_DIR/functions' ]]"
check "utils/ directory exists" "[[ -d '$REPO_DIR/utils' ]]"
check "config1.yaml exists" "[[ -f '$REPO_DIR/config1.yaml' ]]"

# Key scripts
echo -e "\n${BLUE}Key Scripts${NC}"
check "setup_data.py exists" "[[ -f '$REPO_DIR/utils/setup_data.py' ]]"
check "create_data_package.sh exists" "[[ -f '$REPO_DIR/utils/create_data_package.sh' ]]"
check "pipeline_groups_new_mask.py exists" "[[ -f '$REPO_DIR/notebooks/pipeline_groups_new_mask.py' ]]"

# Documentation
echo -e "\n${BLUE}Documentation${NC}"
check "SETUP.md exists" "[[ -f '$REPO_DIR/SETUP.md' ]]"
check "IMPLEMENTATION_SUMMARY.md exists" "[[ -f '$REPO_DIR/IMPLEMENTATION_SUMMARY.md' ]]"
check "OSF_WORKFLOW_CHECKLIST.md exists" "[[ -f '$REPO_DIR/docs/OSF_WORKFLOW_CHECKLIST.md' ]]"
check "OSF_DATA_PACKAGING.md exists" "[[ -f '$REPO_DIR/docs/OSF_DATA_PACKAGING.md' ]]"

# Notebook content
echo -e "\n${BLUE}Notebook Content${NC}"
check "Pre-Step -2 in notebook" "grep -q 'PRE-STEP -2' '$REPO_DIR/notebooks/pipeline_groups_new_mask.py'"
check "Pre-Step -1a in notebook" "grep -q 'PRE-STEP -1a' '$REPO_DIR/notebooks/pipeline_groups_new_mask.py'"
check "Pre-Step -1b in notebook" "grep -q 'PRE-STEP -1b' '$REPO_DIR/notebooks/pipeline_groups_new_mask.py'"
check "Pre-Step -1c in notebook" "grep -q 'PRE-STEP -1c' '$REPO_DIR/notebooks/pipeline_groups_new_mask.py'"

# Platform detection
echo -e "\n${BLUE}Platform Detection${NC}"
check "Platform detection code present" "grep -q 'platform.system()' '$REPO_DIR/notebooks/pipeline_groups_new_mask.py'"
check "macOS path configuration" "grep -q '/Users/glerma/soft/mbrain' '$REPO_DIR/notebooks/pipeline_groups_new_mask.py'"
check "Linux path configuration" "grep -q '/export/home/glerma/glerma/soft/mbrain' '$REPO_DIR/notebooks/pipeline_groups_new_mask.py'"

# Optional: Check for actual data (won't fail if missing)
echo -e "\n${BLUE}Optional: Local Data (Not required for setup)${NC}"
if [[ -d "$REPO_DIR/local" ]]; then
    check "local/ directory found" "[[ -d '$REPO_DIR/local' ]]"
else
    echo -e "${YELLOW}ℹ️${NC}  local/ directory not found (will be auto-downloaded)"
fi

# Summary
echo ""
echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"
total=$((pass + fail))
echo -e "${BLUE}Summary${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}Passed: $pass/$total${NC}"
if [[ $fail -gt 0 ]]; then
    echo -e "${RED}Failed: $fail/$total${NC}"
    exit 1
else
    echo ""
    echo -e "${GREEN}✅ All checks passed! Repository is ready.$NC"
    echo ""
    echo "Next steps:"
    echo ""
    echo "1. Create data package (on machine with data access):"
    echo "   bash utils/create_data_package.sh"
    echo ""
    echo "2. Upload ZIP to OSF.io:"
    echo "   Follow docs/OSF_WORKFLOW_CHECKLIST.md (Phase 2)"
    echo ""
    echo "3. Update OSF URL in setup_data.py:"
    echo "   nano utils/setup_data.py"
    echo "   (change DEFAULT_OSF_URL)"
    echo ""
    echo "4. Test auto-download:"
    echo "   python utils/setup_data.py --check-only"
    echo ""
    exit 0
fi
