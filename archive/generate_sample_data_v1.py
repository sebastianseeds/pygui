#!/usr/bin/env python3
"""
Generate sample data for testing the Electron Flux Visualizer
Creates:
1. A synthetic VTK file with electron flux field
2. A synthetic orbital CSV with LEO satellite trajectory
"""

import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import numpy_to_vtk

def create_sample_vtk_data(filename="sample_electron_flux.vtk"):
    """Create a sample VTK file with synthetic electron flux data"""
    
    print(f"Creating sample VTK data: {filename}")
    
    # Create a 3D grid around Earth (smaller for testing)
    # Grid from -15000 to 15000 km in each direction
    nx, ny, nz = 30, 30, 30  # Reduced size for faster loading
    
    # Create structured grid points
    x = np.linspace(-15000, 15000, nx)
    y = np.linspace(-15000, 15000, ny) 
    z = np.linspace(-15000, 15000, nz)
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Calculate distance from Earth center
    R = np.sqrt(X**2 + Y**2 + Z**2)
    
    # Create synthetic electron flux field
    earth_radius = 6371  # km
    
    # Van Allen belt-like structure
    flux = np.zeros_like(R)
    
    # Inner belt (1.2 - 3 Earth radii)
    inner_belt = (R >= 1.2 * earth_radius) & (R <= 3 * earth_radius)
    flux[inner_belt] = 1e6 * np.exp(-(R[inner_belt] - 2*earth_radius)**2 / (0.5*earth_radius)**2)
    
    # Outer belt (3 - 7 Earth radii) 
    outer_belt = (R >= 3 * earth_radius) & (R <= 7 * earth_radius)
    flux[outer_belt] = 5e5 * np.exp(-(R[outer_belt] - 4*earth_radius)**2 / (earth_radius)**2)
    
    # Add some noise and asymmetry
    np.random.seed(42)  # For reproducible results
    flux += 1e4 * np.random.random(flux.shape)
    flux *= (1 + 0.3 * np.cos(2 * np.arctan2(Y, X)))  # Day/night asymmetry
    
    # Avoid flux inside Earth
    flux[R < earth_radius] = 0
    
    # Create unstructured grid instead of structured grid
    points = vtk.vtkPoints()
    flux_values = []
    
    # Only add points where flux > 0 to reduce file size
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if flux[i,j,k] > 1000:  # Only include significant flux regions
                    points.InsertNextPoint(X[i,j,k], Y[i,j,k], Z[i,j,k])
                    flux_values.append(flux[i,j,k])
    
    print(f"  Reduced to {points.GetNumberOfPoints()} significant flux points")
    
    # Create unstructured grid
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    
    # Create vertex cells (one for each point)
    for i in range(points.GetNumberOfPoints()):
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, i)
        ugrid.InsertNextCell(vertex.GetCellType(), vertex.GetPointIds())
    
    # Add flux data
    flux_array = numpy_to_vtk(np.array(flux_values), deep=True)
    flux_array.SetName("electron_flux")
    ugrid.GetPointData().SetScalars(flux_array)
    
    # Write to file
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(ugrid)
    writer.Write()
    
    print(f"  Points with significant flux: {points.GetNumberOfPoints()}")
    print(f"  Flux range: {min(flux_values):.2e} to {max(flux_values):.2e} particles/cm²/s/sr")
    print(f"  File saved: {filename}")
    
    return ugrid


# def create_sample_vtk_data(filename="sample_electron_flux.vtk"):
#     """Create a sample VTK file with synthetic electron flux data"""
    
#     print(f"Creating sample VTK data: {filename}")
    
#     # Create a 3D grid around Earth
#     # Grid from -15000 to 15000 km in each direction
#     nx, ny, nz = 50, 50, 50
    
#     # Create structured grid
#     x = np.linspace(-15000, 15000, nx)
#     y = np.linspace(-15000, 15000, ny) 
#     z = np.linspace(-15000, 15000, nz)
    
#     X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
#     # Calculate distance from Earth center
#     R = np.sqrt(X**2 + Y**2 + Z**2)
    
#     # Create synthetic electron flux field
#     # Higher flux in radiation belts (~3-7 Earth radii)
#     earth_radius = 6371  # km
    
#     # Van Allen belt-like structure
#     flux = np.zeros_like(R)
    
#     # Inner belt (1.2 - 3 Earth radii)
#     inner_belt = (R >= 1.2 * earth_radius) & (R <= 3 * earth_radius)
#     flux[inner_belt] = 1e6 * np.exp(-(R[inner_belt] - 2*earth_radius)**2 / (0.5*earth_radius)**2)
    
#     # Outer belt (3 - 7 Earth radii) 
#     outer_belt = (R >= 3 * earth_radius) & (R <= 7 * earth_radius)
#     flux[outer_belt] = 5e5 * np.exp(-(R[outer_belt] - 4*earth_radius)**2 / (earth_radius)**2)
    
#     # Add some noise and asymmetry
#     flux += 1e4 * np.random.random(flux.shape)
#     flux *= (1 + 0.3 * np.cos(2 * np.arctan2(Y, X)))  # Day/night asymmetry
    
#     # Avoid flux inside Earth
#     flux[R < earth_radius] = 0
    
#     # Create VTK structured grid
#     points = vtk.vtkPoints()
#     for k in range(nz):
#         for j in range(ny):
#             for i in range(nx):
#                 points.InsertNextPoint(X[i,j,k], Y[i,j,k], Z[i,j,k])
    
#     # Create grid
#     grid = vtk.vtkStructuredGrid()
#     grid.SetDimensions(nx, ny, nz)
#     grid.SetPoints(points)
    
#     # Add flux data
#     flux_vtk = numpy_to_vtk(flux.ravel(), deep=True)
#     flux_vtk.SetName("electron_flux")
#     grid.GetPointData().SetScalars(flux_vtk)
    
#     # Write to file
#     writer = vtk.vtkStructuredGridWriter()
#     writer.SetFileName(filename)
#     writer.SetInputData(grid)
#     writer.Write()
    
#     print(f"  Grid size: {nx}x{ny}x{nz} = {nx*ny*nz} points")
#     print(f"  Flux range: {flux.min():.2e} to {flux.max():.2e} particles/cm²/s/sr")
#     print(f"  File saved: {filename}")

def create_sample_orbital_data(filename="sample_orbit.csv", orbit_type="LEO"):
    """Create sample orbital CSV data"""
    
    print(f"Creating sample orbital data: {filename} ({orbit_type})")
    
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
    
    # Orbital parameters
    inclination_rad = np.radians(inclination)
    
    positions = []
    velocities = []
    
    for t in times:
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
    
    # Save to CSV
    df.to_csv(filename, index=False)
    
    print(f"  Orbit type: {orbit_type}")
    print(f"  Altitude: {altitude} km")
    print(f"  Period: {period} hours")
    print(f"  Inclination: {inclination}°")
    print(f"  Points: {n_points}")
    print(f"  File saved: {filename}")

def main():
    """Generate all sample data files"""
    print("Generating sample data for Electron Flux Visualizer")
    print("=" * 50)
    
    # Create VTK data
    create_sample_vtk_data("sample_electron_flux.vtk")
    print()
    
    # Create different orbital data
    create_sample_orbital_data("sample_leo_orbit.csv", "LEO")
    print()
    create_sample_orbital_data("sample_geo_orbit.csv", "GEO")
    print()
    create_sample_orbital_data("sample_meo_orbit.csv", "MEO")
    print()
    
    print("Sample data generation complete!")
    print()
    print("To test the visualizer:")
    print("1. python flux_visualizer.py")
    print("2. Load VTK Data: sample_electron_flux.vtk")
    print("3. Load Orbital CSV: sample_leo_orbit.csv (or others)")
    print("4. Click Play to start animation")

if __name__ == '__main__':
    main()
