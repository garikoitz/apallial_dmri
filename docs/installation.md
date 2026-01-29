# Installation Guide

## Required Dependencies

### Core Tools

1. **Python 3.11+**
   ```bash
   # Check version
   python --version
   ```

2. **MRtrix3** (for DWI processing)
   ```bash
   # Check installation
   mrconvert --version
   
   # Install on macOS
   brew install mrtrix3
   
   # Or from source: https://www.mrtrix.org/download/
   ```

3. **ANTs** (Advanced Normalization Tools - for coregistration)
   ```bash
   # Check installation
   antsRegistrationSyN.sh --version
   ConvertTransformFile --version
   
   # Install on macOS
   brew install ants
   
   # Or using conda
   conda install -c conda-forge ants
   
   # Or from source: https://github.com/ANTsX/ANTs
   ```

### Python Packages

Install required Python packages:
```bash
cd /Users/glerma/soft/mbrain
pip install -r requirements.txt
```

## Verify Installation

Run the verification script to check all dependencies:

```bash
python -c "
import subprocess
import sys

def check_command(cmd, name):
    try:
        result = subprocess.run([cmd, '--version'], capture_output=True, text=True)
        print(f'✓ {name}: installed')
        return True
    except FileNotFoundError:
        print(f'✗ {name}: NOT FOUND')
        return False

print('Checking dependencies...')
print()

checks = [
    ('mrconvert', 'MRtrix3'),
    ('antsRegistrationSyN.sh', 'ANTs'),
    ('ConvertTransformFile', 'ANTs Tools'),
]

all_ok = True
for cmd, name in checks:
    if not check_command(cmd, name):
        all_ok = False

print()
if all_ok:
    print('✅ All dependencies installed!')
else:
    print('❌ Some dependencies are missing. See above.')
    sys.exit(1)
"
```

## Installation Steps

### 1. Install Homebrew (macOS)

If you don't have Homebrew installed:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### 2. Install MRtrix3

```bash
brew install mrtrix3
```

Verify:
```bash
mrconvert --version
# Should show: mrconvert 3.0.2-48-g63d94827
```

### 3. Install ANTs

```bash
brew install ants
```

Verify:
```bash
antsRegistrationSyN.sh --version
# Should show ANTs version information
```

### 4. Install Python Dependencies

```bash
cd /Users/glerma/soft/mbrain
pip install -r requirements.txt
```

## Test the Pipeline

Run a dry-run test on the first subject:

```bash
python main.py --preprocess --steps coreg --first-subject --dry-run
```

This should show all commands without executing them.

## Troubleshooting

### MRtrix3 Not Found

If MRtrix3 commands are not found, check your PATH:
```bash
echo $PATH
which mrconvert
```

Add to your `~/.zshrc` or `~/.bashrc`:
```bash
export PATH="/opt/homebrew/bin:$PATH"
```

### ANTs Not Found

If ANTs commands are not found:
```bash
# Find ANTs installation
brew list ants

# Add to PATH if needed
export PATH="/opt/homebrew/opt/ants/bin:$PATH"
```

### Python Module Errors

If you get import errors:
```bash
# Ensure you're using the right Python
which python
python --version

# Reinstall packages
pip install --upgrade pip
pip install -r requirements.txt
```

## Optional Tools

### For visualization:
- **ITK-SNAP** (ROI drawing): https://www.itksnap.org/
- **FSLeyes** (MRI viewing): https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes

### For ROI drawing in MRtrix3:
```bash
# MRtrix3 includes mrview with built-in ROI editor
mrview --help
```

## Running the Pipeline

Once all dependencies are installed:

1. **BIDS conversion**:
   ```bash
   python main.py --bids
   ```

2. **Preprocessing** (includes two-stage coregistration):
   ```bash
   python main.py --preprocess
   ```

3. **Tractography**:
   ```bash
   python main.py --tracking
   ```

4. **All steps**:
   ```bash
   python main.py --all
   ```

## Getting Help

```bash
python main.py --help
```
