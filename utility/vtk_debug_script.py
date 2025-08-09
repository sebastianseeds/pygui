#!/usr/bin/env python3
"""
VTK File Diagnostic Script
Analyzes VTK files to debug loading issues
"""

import vtk
import os
from pathlib import Path

def diagnose_vtk_file(filename):
    """Comprehensive VTK file diagnosis"""
    
    print(f"Diagnosing VTK file: {filename}")
    print("=" * 60)
    
    # Check if file exists
    if not os.path.exists(filename):
        print(f"❌ ERROR: File '{filename}' does not exist!")
        return False
    
    file_size = os.path.getsize(filename)
    print(f"File size: {file_size:,} bytes ({file_size/1024/1024:.2f} MB)")
    
    # Read first few lines to check format
    print("\nFile header (first 10 lines):")
    print("-" * 40)
    try:
        with open(filename, 'r') as f:
            for i, line in enumerate(f):
                if i >= 10:
                    break
                print(f"{i+1:2d}: {line.rstrip()}")
    except UnicodeDecodeError:
        print("   (Binary VTK file - cannot display header)")
    
    # Try different VTK readers
    readers = [
        ("UnstructuredGridReader", vtk.vtkUnstructuredGridReader),
        ("StructuredGridReader", vtk.vtkStructuredGridReader),
        ("PolyDataReader", vtk.vtkPolyDataReader),
        ("StructuredPointsReader", vtk.vtkStructuredPointsReader),
        ("RectilinearGridReader", vtk.vtkRectilinearGridReader),
    ]
    
    print(f"\nTesting different VTK readers:")
    print("-" * 40)
    
    successful_readers = []
    
    for reader_name, reader_class in readers:
        try:
            print(f"   Testing {reader_name}...", end=" ")
            
            reader = reader_class()
            reader.SetFileName(filename)
            reader.Update()
            
            output = reader.GetOutput()
            
            if output and output.GetNumberOfPoints() > 0:
                print(f"✅ SUCCESS")
                successful_readers.append((reader_name, reader, output))
                
                # Get details
                n_points = output.GetNumberOfPoints()
                n_cells = output.GetNumberOfCells()
                
                print(f"      Points: {n_points:,}")
                print(f"      Cells: {n_cells:,}")
                
                # Check scalar data
                point_data = output.GetPointData()
                if point_data:
                    n_arrays = point_data.GetNumberOfArrays()
                    print(f"      Point data arrays: {n_arrays}")
                    
                    for i in range(n_arrays):
                        array = point_data.GetArray(i)
                        if array:
                            array_name = array.GetName() or f"Array_{i}"
                            array_size = array.GetNumberOfTuples()
                            array_components = array.GetNumberOfComponents()
                            print(f"        [{i}] {array_name}: {array_size} tuples, {array_components} components")
                            
                            if array_size > 0:
                                val_range = array.GetRange()
                                print(f"            Range: {val_range[0]:.2e} to {val_range[1]:.2e}")
                else:
                    print(f"      No point data found")
                    
            else:
                print(f"No data")
                
        except Exception as e:
            print(f"ERROR: {str(e)}")
    
    if successful_readers:
        print(f"\nFound {len(successful_readers)} working reader(s)")
        
        # Test the conversion process our app uses
        print(f"\nTesting conversion to UnstructuredGrid:")
        print("-" * 40)
        
        for reader_name, reader, output in successful_readers:
            print(f"Converting from {reader_name}...")
            
            try:
                if isinstance(output, vtk.vtkUnstructuredGrid):
                    print(f"      Already UnstructuredGrid - no conversion needed")
                    final_data = output
                    
                elif isinstance(output, vtk.vtkStructuredGrid):
                    print(f"      Converting StructuredGrid...")
                    
                    # Method 1: Direct conversion
                    converter = vtk.vtkStructuredGridGeometryFilter()
                    converter.SetInputData(output)
                    converter.Update()
                    
                    geometry = converter.GetOutput()
                    print(f"       Geometry points: {geometry.GetNumberOfPoints()}")
                    
                    # Convert to unstructured grid
                    append_filter = vtk.vtkAppendFilter()
                    append_filter.AddInputData(geometry)
                    append_filter.Update()
                    
                    final_data = append_filter.GetOutput()
                    print(f"      Converted successfully")
                    
                elif isinstance(output, vtk.vtkPolyData):
                    print(f"      Converting PolyData...")
                    
                    append_filter = vtk.vtkAppendFilter()
                    append_filter.AddInputData(output)
                    append_filter.Update()
                    
                    final_data = append_filter.GetOutput()
                    print(f"      Converted successfully")
                    
                else:
                    print(f"      Unknown data type: {type(output)}")
                    continue
                
                # Check final result
                print(f"      Final points: {final_data.GetNumberOfPoints()}")
                print(f"      Final cells: {final_data.GetNumberOfCells()}")
                
                # Check scalar data preservation
                final_scalars = final_data.GetPointData().GetScalars()
                if final_scalars:
                    print(f"      Scalar data preserved: {final_scalars.GetName()}")
                    print(f"      Range: {final_scalars.GetRange()}")
                else:
                    print(f"      Scalar data lost during conversion!")
                    
                    # Try to find scalar data in other arrays
                    point_data = final_data.GetPointData()
                    if point_data.GetNumberOfArrays() > 0:
                        first_array = point_data.GetArray(0)
                        point_data.SetScalars(first_array)
                        print(f"     Using first array as scalars: {first_array.GetName()}")
                    
            except Exception as e:
                print(f"     Conversion failed: {str(e)}")
                import traceback
                traceback.print_exc()
                
    else:
        print(f"\nNo working readers found!")
        
        # Additional diagnostics
        print(f"\nAdditional diagnostics:")
        print("-" * 40)
        
        # Check if it's a binary file
        try:
            with open(filename, 'rb') as f:
                first_bytes = f.read(100)
                print(f"   First bytes (hex): {first_bytes[:20].hex()}")
                
                # Try to decode as text
                try:
                    text_content = first_bytes.decode('utf-8')
                    print(f"   Appears to be text format")
                except:
                    print(f"   Appears to be binary format")
                    
        except Exception as e:
            print(f"   Cannot read file: {e}")
    
    return len(successful_readers) > 0

def main():
    """Test the generated VTK file"""
    
    filename = "sample_electron_flux.vtk"
    
    success = diagnose_vtk_file(filename)
    
    print(f"\n" + "=" * 60)
    if success:
        print(f"DIAGNOSIS COMPLETE - File appears readable")
        print(f" The issue may be in the GUI application's loading logic")
    else:
        print(f"DIAGNOSIS COMPLETE - File has issues")
        print(f" Need to fix the VTK file generation")
    
    print(f"\nRecommendations:")
    if success:
        print(f"   • Check the exact error message in the GUI")
        print(f"   • Add debug prints to the GUI's load_vtk_data method")
        print(f"   • Verify the GUI is using the correct reader")
    else:
        print(f"   • Regenerate the VTK file with corrected format")
        print(f"   • Try a different VTK file format")
        print(f"   • Check VTK installation")

if __name__ == '__main__':
    main()
