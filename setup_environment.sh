#!/bin/bash

echo "Setting up Python Electron Flux Visualizer Environment"

# Create virtual environment
python3 -m venv flux_env
source flux_env/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install requirements
pip install PyQt6 vtk numpy pandas matplotlib

# Verify installation
python -c "import PyQt6; print('PyQt6 installed successfully')"
python -c "import vtk; print(f'VTK version: {vtk.VTK_VERSION}')"
python -c "import numpy; print(f'NumPy version: {numpy.__version__}')"
python -c "import pandas; print(f'Pandas version: {pandas.__version__}')"

echo ""
echo "Environment setup complete!"
echo "To activate environment: source flux_env/bin/activate"
echo "To run application: python flux_visualizer.py"
echo ""
echo "Sample data format:"
echo "VTK: Any VTK unstructured grid with scalar field data"
echo "CSV: columns should include: time,x,y,z (optionally: vx,vy,vz)"
echo "     time in hours, positions in km, velocities in km/s"
