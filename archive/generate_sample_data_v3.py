#!/usr/bin/env python3
"""
Fixed VTK Data Generator for Electron Flux Visualizer
Creates proper structured grid data that loads correctly in the visualizer.
"""

import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import numpy_to_vtk
import time
import sys

def print_progress_bar(iteration, total, prefix='', suffix='', length=50, fill='█', print_end="\r"):
    """Create terminal progress bar"""
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=print_end)
    if iteration == total: 
        print()

def create_structured_grid_vtk(filename="electron_flux_structured.vts"):
    """Create a proper structured grid VTK file (XML format)"""
    
    print(f"Creating structured grid VTK data: {filename}")
    start_time = time.time()
    
    # Create a 3D grid around Earth
    nx, ny, nz = 40, 40, 40
    total_points = nx * ny * nz
    
    print(f"  Generating {nx}×{ny}×{nz} = {total_points:,} grid points...")
    
    # Create structured grid
    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(nx, ny, nz)
    
    # Create points
    points = vtk.vtkPoints()
    
    # Grid bounds (km from Earth center)
    x_min, x_max = -20000, 20000
    y_min, y_max = -20000, 20000  
    z_min, z_max = -20000, 20000
    
    x = np.linspace(x_min, x_max, nx)
    y = np.linspace(y_min, y_max, ny)
    z = np.linspace(z_min, z_max, nz)
    
    flux_values = []
    earth_radius = 6371  # km
    
    print("  Computing flux field...")
    
    point_count = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                xi, yj, zk = x[i], y[j], z[k]
                points.InsertNextPoint(xi, yj, zk)
                
                # Calculate distance from Earth center
                r = np.sqrt(xi**2 + yj**2 + zk**2)
                
                # Van Allen belt electron flux model
                flux = 0.0
                
                if r > earth_radius:
                    # Inner belt (1.2 - 3 Earth radii)
                    if 1.2 * earth_radius <= r <= 3 * earth_radius:
                        flux += 1e6 * np.exp(-((r - 2*earth_radius)**2) / (0.5*earth_radius)**2)
                    
                    # Outer belt (3 - 7 Earth radii)  
                    if 3 * earth_radius <= r <= 7 * earth_radius:
                        flux += 5e5 * np.exp(-((r - 4*earth_radius)**2) / (earth_radius)**2)
                    
                    # Add some noise and asymmetry
                    phi = np.arctan2(yj, xi)
                    flux *= (1 + 0.3 * np.cos(2 * phi))  # Day/night asymmetry
                    flux += np.random.uniform(0, 1e4)  # Background noise
                
                flux_values.append(max(0, flux))  # Ensure non-negative
                point_count += 1
                
                if point_count % 1000 == 0:
                    print_progress_bar(point_count, total_points, prefix='    Computing', 
                                     suffix=f'({point_count:,}/{total_points:,})')
    
    print_progress_bar(total_points, total_points, prefix='    Computing', 
                      suffix=f'({total_points:,}/{total_points:,})')
    
    sgrid.SetPoints(points)
    
    # Add flux data
    print("  Adding flux data to grid...")
    flux_array = numpy_to_vtk(np.array(flux_values), deep=True)
    flux_array.SetName("electron_flux")
    sgrid.GetPointData().SetScalars(flux_array)
    
    # Write XML structured grid file
    print(f"  Writing XML file: {filename}")
    writer = vtk.vtkXMLStructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(sgrid)
    writer.Write()
    
    elapsed_time = time.time() - start_time
    flux_array_np = np.array(flux_values)
    
    print(f"✓ Structured grid VTK file created!")
    print(f"  Processing time: {elapsed_time:.1f} seconds")
    print(f"  Grid dimensions: {nx} × {ny} × {nz}")
    print(f"  Total points: {total_points:,}")
    print(f"  Flux range: {flux_array_np.min():.2e} to {flux_array_np.max():.2e}")
    print(f"  Non-zero flux points: {np.sum(flux_array_np > 0):,}")
    print(f"  File: {filename}")
    
    return sgrid

def create_unstructured_grid_vtk(filename="electron_flux_unstructured.vtu"):
    """Create an unstructured grid VTK file with only significant flux regions"""
    
    print(f"Creating unstructured grid VTK data: {filename}")
    start_time = time.time()
    
    # Create sampling points in spherical coordinates around Earth
    earth_radius = 6371  # km
    
    # Define radial shells
    r_shells = np.linspace(1.1 * earth_radius, 8 * earth_radius, 20)
    
    # Angular sampling
    n_theta = 30  # polar angle
    n_phi = 60    # azimuthal angle
    
    theta = np.linspace(0, np.pi, n_theta)
    phi = np.linspace(0, 2*np.pi, n_phi)
    
    total_points = len(r_shells) * n_theta * n_phi
    print(f"  Generating {len(r_shells)} × {n_theta} × {n_phi} = {total_points:,} points...")
    
    points = vtk.vtkPoints()
    flux_values = []
    
    point_count = 0
    significant_points = 0
    
    print("  Computing flux in spherical shells...")
    
    for i, r in enumerate(r_shells):
        for j, t in enumerate(theta):
            for k, p in enumerate(phi):
                # Convert to Cartesian
                x = r * np.sin(t) * np.cos(p)
                y = r * np.sin(t) * np.sin(p)
                z = r * np.cos(t)
                
                # Calculate Van Allen belt flux
                flux = 0.0
                
                # Inner belt (1.2 - 3 Earth radii)
                if 1.2 * earth_radius <= r <= 3 * earth_radius:
                    flux += 1e6 * np.exp(-((r - 2*earth_radius)**2) / (0.5*earth_radius)**2)
                
                # Outer belt (3 - 7 Earth radii)
                if 3 * earth_radius <= r <= 7 * earth_radius:
                    flux += 5e5 * np.exp(-((r - 4*earth_radius)**2) / (earth_radius)**2)
                
                # Add asymmetry and noise
                if flux > 0:
                    flux *= (1 + 0.3 * np.cos(2 * p))  # Day/night
                    flux *= (1 + 0.1 * np.sin(3 * t))  # Latitudinal variation
                    flux += np.random.uniform(0, 1e4)
                
                # Only keep points with significant flux
                if flux > 1000:  # Threshold for significance
                    points.InsertNextPoint(x, y, z)
                    flux_values.append(flux)
                    significant_points += 1
                
                point_count += 1
                if point_count % 5000 == 0:
                    print_progress_bar(point_count, total_points, prefix='    Sampling',
                                     suffix=f'({significant_points:,} significant)')
    
    print_progress_bar(total_points, total_points, prefix='    Sampling',
                      suffix=f'({significant_points:,} significant)')
    
    print(f"  Creating unstructured grid with {significant_points:,} points...")
    
    # Create unstructured grid
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    
    # Create vertex cells
    for i in range(points.GetNumberOfPoints()):
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, i)
        ugrid.InsertNextCell(vertex.GetCellType(), vertex.GetPointIds())
        
        if (i + 1) % 1000 == 0:
            print_progress_bar(i + 1, points.GetNumberOfPoints(), prefix='    Cells')
    
    print_progress_bar(points.GetNumberOfPoints(), points.GetNumberOfPoints(), prefix='    Cells')
    
    # Add flux data
    print("  Adding flux data...")
    flux_array = numpy_to_vtk(np.array(flux_values), deep=True)
    flux_array.SetName("electron_flux")
    ugrid.GetPointData().SetScalars(flux_array)
    
    # Write XML unstructured grid file
    print(f"  Writing XML file: {filename}")
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(ugrid)
    writer.Write()
    
    elapsed_time = time.time() - start_time
    flux_array_np = np.array(flux_values)
    
    print(f"✓ Unstructured grid VTK file created!")
    print(f"  Processing time: {elapsed_time:.1f} seconds")
    print(f"  Significant points: {significant_points:,}")
    print(f"  Flux range: {flux_array_np.min():.2e} to {flux_array_np.max():.2e}")
    print(f"  File: {filename}")
    
    return ugrid

def create_legacy_vtk(filename="electron_flux_legacy.vtk"):
    """Create a legacy format VTK file for compatibility"""
    
    print(f"Creating legacy VTK file: {filename}")
    start_time = time.time()
    
    # Simple structured points (image data)
    nx, ny, nz = 50, 50, 50
    
    # Create image data
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(nx, ny, nz)
    imageData.SetSpacing(800, 800, 800)  # 800 km spacing
    imageData.SetOrigin(-20000, -20000, -20000)  # Centered on Earth
    
    # Generate flux data
    print("  Generating flux field...")
    flux_values = []
    earth_radius = 6371
    
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # Get world coordinates
                x = -20000 + i * 800
                y = -20000 + j * 800  
                z = -20000 + k * 800
                
                r = np.sqrt(x**2 + y**2 + z**2)
                
                flux = 0.0
                if r > earth_radius:
                    # Simple belt model
                    if 1.5 * earth_radius <= r <= 6 * earth_radius:
                        flux = 1e6 * np.exp(-((r - 3*earth_radius)**2) / (earth_radius)**2)
                        phi = np.arctan2(y, x)
                        flux *= (1 + 0.5 * np.cos(2 * phi))
                
                flux_values.append(max(0, flux))
    
    # Add scalar data
    flux_array = numpy_to_vtk(np.array(flux_values), deep=True)
    flux_array.SetName("electron_flux")
    imageData.GetPointData().SetScalars(flux_array)
    
    # Write legacy VTK file
    print(f"  Writing legacy VTK file: {filename}")
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(filename)
    writer.SetInputData(imageData)
    writer.Write()
    
    elapsed_time = time.time() - start_time
    flux_array_np = np.array(flux_values)
    
    print(f"✓ Legacy VTK file created!")
    print(f"  Processing time: {elapsed_time:.1f} seconds")
    print(f"  Dimensions: {nx} × {ny} × {nz}")
    print(f"  Spacing: 800 km")
    print(f"  Flux range: {flux_array_np.min():.2e} to {flux_array_np.max():.2e}")
    print(f"  File: {filename}")

def create_sample_orbital_data(filename="sample_orbit.csv", orbit_type="LEO"):
    """Create sample orbital CSV data - same as before but with fixes"""
    
    print(f"Creating orbital data: {filename} ({orbit_type})")
    start_time = time.time()
    
    if orbit_type == "LEO":
        altitude = 400  # km
        period = 1.5    # hours
        inclination = 51.6  # degrees
    elif orbit_type == "GEO":
        altitude = 35786  # km
        period = 24.0   # hours  
        inclination = 0.0  # degrees
    else:  # MEO
        altitude = 20200  # km
        period = 12.0   # hours
        inclination = 55.0  # degrees
    
    earth_radius = 6371  # km
    orbital_radius = earth_radius + altitude
    
    # More points for smoother animation
    n_points = 300
    times = np.linspace(0, period, n_points)
    
    print(f"  Computing {n_points} orbital positions...")
    
    inclination_rad = np.radians(inclination)
    
    positions = []
    velocities = []
    
    for idx, t in enumerate(times):
        # Mean anomaly
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
        
        # Velocity
        orbital_speed = 2 * np.pi * orbital_radius / (period * 3600)  # km/s
        
        vx_orb = -orbital_speed * np.sin(M)
        vy_orb = orbital_speed * np.cos(M)
        vz_orb = 0
        
        vx = vx_orb
        vy = vy_orb * np.cos(inclination_rad) - vz_orb * np.sin(inclination_rad)
        vz = vy_orb * np.sin(inclination_rad) + vz_orb * np.cos(inclination_rad)
        
        velocities.append([vx, vy, vz])
        
        if (idx + 1) % 50 == 0:
            print_progress_bar(idx + 1, n_points, prefix='    Computing')
    
    print_progress_bar(n_points, n_points, prefix='    Computing')
    
    positions = np.array(positions)
    velocities = np.array(velocities)
    
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
    
    df.to_csv(filename, index=False)
    
    elapsed_time = time.time() - start_time
    print(f"  ✓ Orbital data created! ({elapsed_time:.2f}s)")
    print(f"    Points: {n_points}, Period: {period}h, Alt: {altitude:,}km")

def main():
    """Generate multiple VTK formats for testing"""
    print("Electron Flux Visualizer - Enhanced Data Generator")
    print("=" * 60)
    
    overall_start = time.time()
    
    print("\nGenerating VTK files in different formats...")
    print("-" * 50)
    
    # Create XML structured grid (recommended)
    print("\n1. XML Structured Grid (.vts)")
    create_structured_grid_vtk("flux_field_structured.vts")
    
    # Create XML unstructured grid  
    print("\n2. XML Unstructured Grid (.vtu)")
    create_unstructured_grid_vtk("flux_field_unstructured.vtu")
    
    # Create legacy VTK
    print("\n3. Legacy VTK (.vtk)")
    create_legacy_vtk("flux_field_legacy.vtk")
    
    # Create orbital data
    print("\n" + "-" * 50)
    print("Generating orbital data...")
    create_sample_orbital_data("orbit_leo.csv", "LEO")
    
    total_time = time.time() - overall_start
    
    print("\n" + "=" * 60)
    print("✓ All sample data generated successfully!")
    print(f"  Total time: {total_time:.1f} seconds")
    print("\nGenerated files:")
    print("  • flux_field_structured.vts   - XML structured grid (RECOMMENDED)")
    print("  • flux_field_unstructured.vtu - XML unstructured grid")  
    print("  • flux_field_legacy.vtk       - Legacy VTK format")
    print("  • orbit_leo.csv               - LEO satellite trajectory")
    
    print("\nTo test:")
    print("  1. python flux_visualizer.py")
    print("  2. Load VTK Data → flux_field_structured.vts")
    print("  3. Load Orbital CSV → orbit_leo.csv")
    print("  4. Click Play to animate")

if __name__ == '__main__':
    main()
