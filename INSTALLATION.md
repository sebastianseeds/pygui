# Installation Guide

## Quick Setup (Recommended)

Use the provided shell script for automatic setup:

```bash
chmod +x setup_environment.sh
./setup_environment.sh
```

This will create a virtual environment, install all dependencies, and verify the installation.

## Manual Installation

### Step 1: Create Virtual Environment
```bash
python3 -m venv flux_env
source flux_env/bin/activate  # Linux/macOS
# flux_env\Scripts\activate     # Windows
```

### Step 2: Install Dependencies
```bash
pip install --upgrade pip
pip install numpy pandas matplotlib PyQt6 vtk
```

### Step 3: Verify Installation
```bash
python -c "import vtk; print(f'VTK version: {vtk.VTK_VERSION}')"
python -c "import PyQt6; print('PyQt6 OK')"
```

## VTK Installation Troubleshooting

VTK can be challenging to install. Try these methods in order:

### Method 1: pip (standard)
```bash
pip install vtk
```

### Method 2: conda (most reliable)
```bash
conda install -c conda-forge vtk
```

### Method 3: System packages

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install python3-vtk9 python3-vtk9-dev
```

**macOS (Homebrew):**
```bash
brew install vtk
pip install vtk
```

**Fedora:**
```bash
sudo dnf install vtk-devel python3-vtk
```

**CentOS/RHEL:**
```bash
sudo yum install vtk-devel python3-vtk
```

## Common Issues & Solutions

### "No module named 'vtk'"
VTK installation failed. Try the conda method:
```bash
conda install -c conda-forge vtk
```

### "Qt platform plugin could not be initialized"
Set the Qt platform:
```bash
export QT_QPA_PLATFORM=xcb     # Linux
export QT_QPA_PLATFORM=cocoa   # macOS
```

### "ImportError: libGL.so.1" (Linux)
Install OpenGL libraries:
```bash
sudo apt-get install libgl1-mesa-glx libglu1-mesa
```

### Permission denied on shell script
Make it executable:
```bash
chmod +x setup_environment.sh
```

## Requirements

- **Python 3.8+**
- **Operating System:** Linux, macOS, or Windows
- **Memory:** 4GB+ RAM recommended for large datasets
- **Graphics:** OpenGL-compatible graphics card

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| PyQt6 | ≥6.0.0 | GUI framework |
| VTK | ≥9.0.0 | 3D visualization |
| NumPy | ≥1.20.0 | Numerical computing |
| Pandas | ≥1.3.0 | Data handling |
| Matplotlib | ≥3.5.0 | 2D plotting |

## Running the Application

After successful installation:

```bash
# Activate environment
source flux_env/bin/activate

# Run the application
python flux_visualizer.py
```

## Data Format Requirements

### VTK Files
- **Formats:** .vtk, .vtu, .vts, .vtp, .vti
- **Content:** Scalar field data representing electron flux
- **Structure:** Structured or unstructured grids supported

### CSV Files (Orbital Data)
- **Required columns:** `time`, `x`, `y`, `z`
- **Optional columns:** `vx`, `vy`, `vz`
- **Units:** 
  - Time: hours
  - Positions: kilometers
  - Velocities: km/s

**Example CSV format:**
```csv
time,x,y,z,vx,vy,vz
0.0,-6000,2000,1000,5.2,-1.1,0.3
0.1,-5950,2100,1050,5.1,-1.0,0.4
```

## Getting Help

If you encounter issues:

1. Check this troubleshooting guide
2. Verify your Python version: `python3 --version`
3. Try the conda installation method for VTK
4. Check VTK documentation: https://vtk.org/download/

## Development Setup

For development work:

```bash
# Clone with development dependencies
pip install -e .
pip install pytest black flake8  # Optional dev tools
```