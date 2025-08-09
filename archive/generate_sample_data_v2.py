#!/usr/bin/env python3
"""
Generate sample data for testing the Electron Flux Visualizer
Creates:
1. A synthetic VTK file with electron flux field
2. A synthetic orbital CSV with LEO satellite trajectory

Enhanced with progress bars and timing information.
"""

import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import numpy_to_vtk
import time
import sys

def print_progress_bar(iteration, total, prefix='', suffix='', length=50, fill='█', print_end="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        print_end   - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=print_end)
    # Print New Line on Complete
    if iteration == total: 
        print()

def create_sample_vtk_data(filename="sample_electron_flux.vtk"):
    """Create a sample VTK file with synthetic electron flux data"""
    
    print(f"Creating sample VTK data: {filename}")
    start_time = time.time()
    
    # Create a 3D grid around Earth (smaller for testing)
    nx, ny, nz = 30, 30, 30  # Reduced size for faster loading
    total_points = nx * ny * nz
    
    print(f"  Generating {nx}×{ny}×{nz} = {total_points:,} grid points...")
    
    # Create structured grid points
    x = np.linspace(-15000, 15000, nx)
    y = np.linspace(-15000, 15000, ny) 
    z = np.linspace(-15000, 15000, nz)
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    print("  Computing distances from Earth center...")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    
    print("  Generating Van Allen belt electron flux field...")
    earth_radius = 6371  # km
    
    # Van Allen belt-like structure
    flux = np.zeros_like(R)
    
    # Progress tracking for flux calculation
    print("  Creating inner radiation belt...")
    inner_belt = (R >= 1.2 * earth_radius) & (R <= 3 * earth_radius)
    flux[inner_belt] = 1e6 * np.exp(-(R[inner_belt] - 2*earth_radius)**2 / (0.5*earth_radius)**2)
    
    print("  Creating outer radiation belt...")
    outer_belt = (R >= 3 * earth_radius) & (R <= 7 * earth_radius)
    flux[outer_belt] = 5e5 * np.exp(-(R[outer_belt] - 4*earth_radius)**2 / (earth_radius)**2)
    
    print("  Adding noise and day/night asymmetry...")
    np.random.seed(42)  # For reproducible results
    flux += 1e4 * np.random.random(flux.shape)
    flux *= (1 + 0.3 * np.cos(2 * np.arctan2(Y, X)))  # Day/night asymmetry
    
    # Avoid flux inside Earth
    flux[R < earth_radius] = 0
    
    print("  Building VTK data structure...")
    
    # Create unstructured grid
    points = vtk.vtkPoints()
    flux_values = []
    
    # Count significant points first for progress bar
    significant_mask = flux > 1000
    significant_count = np.sum(significant_mask)
    
    print(f"  Found {significant_count:,} points with significant flux")
    print("  Adding points to VTK structure...")
    
    # Add points with progress tracking
    added_points = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if significant_mask[i,j,k]:
                    points.InsertNextPoint(X[i,j,k], Y[i,j,k], Z[i,j,k])
                    flux_values.append(flux[i,j,k])
                    added_points += 1
                    
                    # Update progress every 100 points
                    if added_points % 100 == 0 or added_points == significant_count:
                        print_progress_bar(
                            added_points, significant_count,
                            prefix='    Points', 
                            suffix=f'({added_points:,}/{significant_count:,})'
                        )
    
    print(f"  Creating unstructured grid with {points.GetNumberOfPoints()} points...")
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    
    # Create vertex cells with progress
    print("  Creating VTK cells...")
    num_points = points.GetNumberOfPoints()
    for i in range(num_points):
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, i)
        ugrid.InsertNextCell(vertex.GetCellType(), vertex.GetPointIds())
        
        # Progress every 1000 cells
        if (i + 1) % 1000 == 0 or i == num_points - 1:
            print_progress_bar(
                i + 1, num_points,
                prefix='    Cells',
                suffix=f'({i+1:,}/{num_points:,})'
            )
    
    print("  Adding flux data to VTK structure...")
    flux_array = numpy_to_vtk(np.array(flux_values), deep=True)
    flux_array.SetName("electron_flux")
    ugrid.GetPointData().SetScalars(flux_array)
    
    print(f"  Writing to file: {filename}")
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(ugrid)
    writer.Write()
    
    elapsed_time = time.time() - start_time
    file_size = vtk.vtkFileOutputWindow.GetInstance()
    
    print(f"✓ VTK file created successfully!")
    print(f"  Processing time: {elapsed_time:.1f} seconds")
    print(f"  Final points: {points.GetNumberOfPoints():,}")
    print(f"  Flux range: {min(flux_values):.2e} to {max(flux_values):.2e} particles/cm²/s/sr")
    print(f"  File: {filename}")
    
    return ugrid

def create_sample_orbital_data(filename="sample_orbit.csv", orbit_type="LEO"):
    """Create sample orbital CSV data with progress tracking"""
    
    print(f"Creating sample orbital data: {filename} ({orbit_type})")
    start_time = time.time()
    
    if orbit_type == "LEO":
        # Low Earth Orbit (ISS-like)
        altitude = 400  # km above Earth surface
        period = 1.5    # hours
        inclination = 51.6  # degrees
        
    elif orbit_type == "GEO":
        # Geostationary orbit
        altitude = 35786  # km
        period = 24.0   # hours  
        inclination = 0.0  # degrees
        
    else:  # MEO
        # Medium Earth Orbit (GPS-like)
        altitude = 20200  # km
        period = 12.0   # hours
        inclination = 55.0  # degrees
    
    earth_radius = 6371  # km
    orbital_radius = earth_radius + altitude
    
    # Time array for one complete orbit
    n_points = 200
    times = np.linspace(0, period, n_points)
    
    print(f"  Computing {n_points} orbital positions...")
    
    # Orbital parameters
    inclination_rad = np.radians(inclination)
    
    positions = []
    velocities = []
    
    for idx, t in enumerate(times):
        # Mean anomaly (simplified circular orbit)
        M = 2 * np.pi * t / period
        
        # Position in orbital plane
        x_orb = orbital_radius * np.cos(M)
        y_orb = orbital_radius * np.sin(M)
        z_orb = 0
        
        # Rotate by inclination
        x = x_orb
        y = y_orb * np.cos(inclination_rad) - z_orb * np.sin(inclination_rad)
        z = y_orb * np.sin(inclination_rad) + z_orb * np.cos(inclination_rad)
        
        positions.append([x, y, z])
        
        # Velocity (tangential to orbit)
        orbital_speed = 2 * np.pi * orbital_radius / (period * 3600)  # km/s
        
        vx_orb = -orbital_speed * np.sin(M)
        vy_orb = orbital_speed * np.cos(M)
        vz_orb = 0
        
        # Rotate velocity by inclination  
        vx = vx_orb
        vy = vy_orb * np.cos(inclination_rad) - vz_orb * np.sin(inclination_rad)
        vz = vy_orb * np.sin(inclination_rad) + vz_orb * np.cos(inclination_rad)
        
        velocities.append([vx, vy, vz])
        
        # Progress bar
        print_progress_bar(
            idx + 1, n_points,
            prefix='    Computing',
            suffix=f'({idx+1}/{n_points})'
        )
    
    positions = np.array(positions)
    velocities = np.array(velocities)
    
    print("  Creating DataFrame and saving to CSV...")
    
    # Create DataFrame
    df = pd.DataFrame({
        'time': times,
        'x': positions[:, 0],
        'y': positions[:, 1], 
        'z': positions[:, 2],
        'vx': velocities[:, 0],
        'vy': velocities[:, 1],
        'vz': velocities[:, 2]
    })
    
    # Save to CSV
    df.to_csv(filename, index=False)
    
    elapsed_time = time.time() - start_time
    
    print(f"  Orbital data created successfully!")
    print(f"  Processing time: {elapsed_time:.2f} seconds")
    print(f"  Orbit type: {orbit_type}")
    print(f"  Altitude: {altitude:,} km")
    print(f"  Period: {period} hours")
    print(f"  Inclination: {inclination}°")
    print(f"  Points: {n_points}")
    print(f"  File: {filename}")

def main():
    """Generate all sample data files with progress tracking"""
    print("Electron Flux Visualizer - Sample Data Generator")
    print("=" * 60)
    
    overall_start = time.time()
    
    # Create VTK data
    print("\nSTEP 1: Creating VTK flux field data")
    print("-" * 40)
    create_sample_vtk_data("sample_electron_flux.vtk")
    
    # Create different orbital data
    print("\nSTEP 2: Creating orbital trajectory data")
    print("-" * 40)
    
    create_sample_orbital_data("sample_leo_orbit.csv", "LEO")
    print()
    create_sample_orbital_data("sample_geo_orbit.csv", "GEO")
    print()
    create_sample_orbital_data("sample_meo_orbit.csv", "MEO")
    
    total_time = time.time() - overall_start
    
    print("\n" + "=" * 60)
    print("  Sample data generation complete!")
    print(f" ️  Total processing time: {total_time:.1f} seconds")
    print("\n  Generated files:")
    print("   • sample_electron_flux.vtk - Van Allen belt electron flux field")
    print("   • sample_leo_orbit.csv - ISS-like Low Earth Orbit")
    print("   • sample_geo_orbit.csv - Geostationary Orbit") 
    print("   • sample_meo_orbit.csv - GPS-like Medium Earth Orbit")
    
    print("\n  To test the visualizer:")
    print("   1. python flux_visualizer.py")
    print("   2. Load VTK Data → sample_electron_flux.vtk")
    print("   3. Load Orbital CSV → sample_leo_orbit.csv")
    print("   4. Click 'Show ...' buttons to open plot windows")
    print("   5. Click 'Play' to start animation")

if __name__ == '__main__':
    main()
