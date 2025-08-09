#!/usr/bin/env python3
"""
Electron Flux Orbital Visualizer
A scientific visualization tool for analyzing electron flux data along orbital paths.
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QPushButton, QSlider, QLabel, QDoubleSpinBox, QFileDialog, 
    QMessageBox, QProgressBar, QGroupBox, QFormLayout, QSpinBox,
    QSplitter, QFrame
)
from PyQt6.QtCore import QTimer, Qt, pyqtSignal
from PyQt6.QtGui import QFont

# VTK-Qt integration
try:
    from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
except ImportError:
    from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

class OrbitalPoint:
    """Data structure to hold orbital position data"""
    def __init__(self, time, x, y, z, vx=0, vy=0, vz=0):
        self.time = time  # hours
        self.x, self.y, self.z = x, y, z  # position in km
        self.vx, self.vy, self.vz = vx, vy, vz  # velocity in km/s
        self.phi = np.arctan2(y, x)  # azimuthal angle
        if self.phi < 0:
            self.phi += 2 * np.pi

class FluxAnalyzer:
    """Handles flux analysis calculations"""
    
    def __init__(self):
        self.vtk_data = None
        self.orbital_path = []
        self.cross_section_radius = 1.0  # meters
        
    def set_vtk_data(self, vtk_data):
        """Set the VTK dataset containing flux field"""
        self.vtk_data = vtk_data
        
    def set_orbital_data(self, orbital_points):
        """Set the orbital path data"""
        self.orbital_path = orbital_points
        
    def set_cross_section(self, radius_meters):
        """Set object cross-sectional radius"""
        self.cross_section_radius = radius_meters
        
    def analyze_flux_at_point(self, orbital_point):
        """Analyze flux at a specific orbital point"""
        if not self.vtk_data:
            return 0.0
            
        # Create a probe to sample the field at this point
        probe = vtk.vtkProbeFilter()
        
        # Create a point to probe
        points = vtk.vtkPoints()
        points.InsertNextPoint(orbital_point.x, orbital_point.y, orbital_point.z)
        
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        
        probe.SetInputData(polydata)
        probe.SetSourceData(self.vtk_data)
        probe.Update()
        
        # Get the flux value
        result = probe.GetOutput()
        if result.GetNumberOfPoints() > 0:
            scalar_array = result.GetPointData().GetScalars()
            if scalar_array:
                flux_value = scalar_array.GetValue(0)
                # Calculate integrated flux through cross-sectional area
                area = np.pi * (self.cross_section_radius ** 2)
                return flux_value * area
                
        return 0.0

class SlicePlotWindow(QMainWindow):
    """Window for displaying orbital slice visualization"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Orbital Slice Visualization")
        self.setGeometry(100, 100, 800, 600)
        
        # Create VTK widget
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        self.setCentralWidget(self.vtk_widget)
        
        # VTK pipeline
        self.renderer = vtk.vtkRenderer()
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.renderer.SetBackground(0.1, 0.1, 0.2)  # Dark blue
        
        # Slice plane and cutter
        self.slice_plane = vtk.vtkPlane()
        self.cutter = vtk.vtkCutter()
        self.cutter.SetCutFunction(self.slice_plane)
        
        # Slice actor
        self.slice_mapper = vtk.vtkDataSetMapper()
        self.slice_mapper.SetInputConnection(self.cutter.GetOutputPort())
        
        self.slice_actor = vtk.vtkActor()
        self.slice_actor.SetMapper(self.slice_mapper)
        self.renderer.AddActor(self.slice_actor)
        
        # Object representation (small sphere)
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(50.0)  # 50 km for visibility
        
        object_mapper = vtk.vtkPolyDataMapper()
        object_mapper.SetInputConnection(sphere.GetOutputPort())
        
        self.object_actor = vtk.vtkActor()
        self.object_actor.SetMapper(object_mapper)
        self.object_actor.GetProperty().SetColor(1.0, 0.0, 0.0)  # Red
        self.renderer.AddActor(self.object_actor)
        
        # Scalar bar
        self.scalar_bar = vtk.vtkScalarBarActor()
        self.scalar_bar.SetLookupTable(self.slice_mapper.GetLookupTable())
        self.scalar_bar.SetTitle("Electron Flux")
        self.renderer.AddActor2D(self.scalar_bar)
        
    def update_slice(self, phi_angle, vtk_data):
        """Update the slice at constant phi angle"""
        if vtk_data is None:
            return
            
        # Set slice plane for constant phi
        nx = -np.sin(phi_angle)
        ny = np.cos(phi_angle) 
        nz = 0.0
        
        self.slice_plane.SetNormal(nx, ny, nz)
        self.slice_plane.SetOrigin(0.0, 0.0, 0.0)
        
        self.cutter.SetInputData(vtk_data)
        self.cutter.Update()
        
        self.vtk_widget.GetRenderWindow().Render()
        
    def set_object_position(self, x, y, z):
        """Update object position"""
        self.object_actor.SetPosition(x, y, z)
        self.vtk_widget.GetRenderWindow().Render()

class SpectrumPlotWindow(QMainWindow):
    """Window for energy spectrum plotting"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Energy Spectrum")
        self.setGeometry(200, 200, 600, 400)
        
        # For now, just a placeholder - we'll add matplotlib later
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        layout = QVBoxLayout(central_widget)
        self.status_label = QLabel("Energy spectrum visualization\n(Matplotlib integration coming soon)")
        self.status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.status_label)
        
    def update_spectrum(self, energy_bins, spectrum_data):
        """Update the energy spectrum plot"""
        # Placeholder for matplotlib spectrum plotting
        self.status_label.setText(f"Spectrum updated with {len(spectrum_data)} energy bins")

class FluxTimePlotWindow(QMainWindow):
    """Window for flux vs time plotting"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Flux vs Time")
        self.setGeometry(300, 300, 600, 400)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        layout = QVBoxLayout(central_widget)
        self.status_label = QLabel("Flux vs time visualization\n(Matplotlib integration coming soon)")
        self.status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.status_label)
        
        self.flux_data = []
        self.time_data = []
        
    def add_flux_point(self, time, flux):
        """Add a new flux measurement point"""
        self.time_data.append(time)
        self.flux_data.append(flux)
        self.status_label.setText(f"Flux data: {len(self.flux_data)} points\n"
                                f"Latest: {flux:.2e} particles/s at t={time:.2f}h")
        
    def clear_data(self):
        """Clear all flux data"""
        self.flux_data.clear()
        self.time_data.clear()
        self.status_label.setText("Flux vs time visualization\n(Data cleared)")

class ElectronFluxVisualizer(QMainWindow):
    """Main application window"""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Electron Flux Orbital Visualizer")
        self.setGeometry(100, 100, 1200, 800)
        
        # Data
        self.vtk_data = None
        self.orbital_path = []
        self.current_time_index = 0
        self.is_playing = False
        
        # Analysis
        self.flux_analyzer = FluxAnalyzer()
        
        # Plot windows
        self.slice_window = None
        self.spectrum_window = None
        self.flux_time_window = None
        
        # Animation timer
        self.animation_timer = QTimer()
        self.animation_timer.timeout.connect(self.animation_step)
        
        self.setup_ui()
        self.setup_vtk()
        
    def setup_ui(self):
        """Setup the user interface"""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QHBoxLayout(central_widget)
        
        # Create splitter
        splitter = QSplitter(Qt.Orientation.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left panel - VTK 3D visualization
        vtk_frame = QFrame()
        vtk_layout = QVBoxLayout(vtk_frame)
        
        # VTK widget
        self.vtk_widget = QVTKRenderWindowInteractor(vtk_frame)
        vtk_layout.addWidget(self.vtk_widget)
        
        splitter.addWidget(vtk_frame)
        
        # Right panel - controls
        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)
        
        # File loading group
        file_group = QGroupBox("Data Loading")
        file_layout = QFormLayout(file_group)
        
        self.load_vtk_button = QPushButton("Load VTK Data")
        self.load_orbit_button = QPushButton("Load Orbital CSV")
        
        file_layout.addRow(self.load_vtk_button)
        file_layout.addRow(self.load_orbit_button)
        
        control_layout.addWidget(file_group)
        
        # Animation controls group
        anim_group = QGroupBox("Animation Controls")
        anim_layout = QVBoxLayout(anim_group)
        
        # Play/pause/stop buttons
        button_layout = QHBoxLayout()
        self.play_button = QPushButton("Play")
        self.pause_button = QPushButton("Pause")
        self.stop_button = QPushButton("Stop")
        
        button_layout.addWidget(self.play_button)
        button_layout.addWidget(self.pause_button)
        button_layout.addWidget(self.stop_button)
        anim_layout.addLayout(button_layout)
        
        # Time slider
        self.time_slider = QSlider(Qt.Orientation.Horizontal)
        self.time_slider.setMinimum(0)
        self.time_slider.setMaximum(100)
        anim_layout.addWidget(QLabel("Time:"))
        anim_layout.addWidget(self.time_slider)
        
        self.time_label = QLabel("00:00:00")
        anim_layout.addWidget(self.time_label)
        
        # Speed control
        speed_layout = QFormLayout()
        self.speed_spinbox = QSpinBox()
        self.speed_spinbox.setRange(10, 1000)
        self.speed_spinbox.setValue(100)
        self.speed_spinbox.setSuffix(" ms")
        speed_layout.addRow("Animation Speed:", self.speed_spinbox)
        anim_layout.addLayout(speed_layout)
        
        control_layout.addWidget(anim_group)
        
        # Analysis parameters group
        analysis_group = QGroupBox("Analysis Parameters")
        analysis_layout = QFormLayout(analysis_group)
        
        self.cross_section_spinbox = QDoubleSpinBox()
        self.cross_section_spinbox.setRange(0.1, 100.0)
        self.cross_section_spinbox.setValue(1.0)
        self.cross_section_spinbox.setSingleStep(0.1)
        self.cross_section_spinbox.setSuffix(" m")
        
        analysis_layout.addRow("Cross Section Radius:", self.cross_section_spinbox)
        control_layout.addWidget(analysis_group)
        
        # Plot windows group
        plots_group = QGroupBox("Plot Windows")
        plots_layout = QVBoxLayout(plots_group)
        
        self.show_slice_button = QPushButton("Show Slice Plot")
        self.show_spectrum_button = QPushButton("Show Spectrum Plot")
        self.show_flux_time_button = QPushButton("Show Flux vs Time")
        
        plots_layout.addWidget(self.show_slice_button)
        plots_layout.addWidget(self.show_spectrum_button)
        plots_layout.addWidget(self.show_flux_time_button)
        
        control_layout.addWidget(plots_group)
        
        # Status
        self.status_label = QLabel("Ready - Load VTK data and orbital CSV to begin")
        self.status_label.setWordWrap(True)
        control_layout.addWidget(self.status_label)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        control_layout.addWidget(self.progress_bar)
        
        control_layout.addStretch()
        
        control_panel.setMaximumWidth(350)
        splitter.addWidget(control_panel)
        
        # Set splitter proportions
        splitter.setSizes([800, 350])
        
        self.connect_signals()
        
    def setup_vtk(self):
        """Setup VTK 3D visualization"""
        # Renderer
        self.renderer = vtk.vtkRenderer()
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.renderer.SetBackground(0.05, 0.05, 0.15)  # Dark blue
        
        # Camera setup for Earth-centered view
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(20000, 20000, 10000)  # km
        camera.SetFocalPoint(0, 0, 0)
        camera.SetViewUp(0, 0, 1)
        
        # Create Earth representation (optional)
        self.create_earth_representation()
        
    def create_earth_representation(self):
        """Create a simple Earth sphere for reference"""
        earth_sphere = vtk.vtkSphereSource()
        earth_sphere.SetRadius(6371.0)  # Earth radius in km
        earth_sphere.SetThetaResolution(50)
        earth_sphere.SetPhiResolution(50)
        
        earth_mapper = vtk.vtkPolyDataMapper()
        earth_mapper.SetInputConnection(earth_sphere.GetOutputPort())
        
        self.earth_actor = vtk.vtkActor()
        self.earth_actor.SetMapper(earth_mapper)
        self.earth_actor.GetProperty().SetColor(0.3, 0.3, 0.8)  # Blue-ish
        self.earth_actor.GetProperty().SetOpacity(0.3)
        
        self.renderer.AddActor(self.earth_actor)
        
    def connect_signals(self):
        """Connect UI signals to slots"""
        self.load_vtk_button.clicked.connect(self.load_vtk_data)
        self.load_orbit_button.clicked.connect(self.load_orbital_data)
        
        self.play_button.clicked.connect(self.start_animation)
        self.pause_button.clicked.connect(self.pause_animation)
        self.stop_button.clicked.connect(self.stop_animation)
        
        self.time_slider.valueChanged.connect(self.jump_to_time)
        self.speed_spinbox.valueChanged.connect(self.set_animation_speed)
        self.cross_section_spinbox.valueChanged.connect(self.update_cross_section)
        
        self.show_slice_button.clicked.connect(self.show_slice_window)
        self.show_spectrum_button.clicked.connect(self.show_spectrum_window)
        self.show_flux_time_button.clicked.connect(self.show_flux_time_window)

    def load_vtk_data(self):
        """Load VTK data file with improved error handling"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load VTK Data", "", "VTK Files (*.vtk *.vtu *.vtp *.vts);;All Files (*)"
        )
        
        if not file_path:
            return
            
        try:
            # Determine file type and use appropriate reader
            file_ext = Path(file_path).suffix.lower()
            
            if file_ext == '.vtk':
                # Legacy VTK format - need to detect type
                reader = self.create_legacy_vtk_reader(file_path)
            elif file_ext == '.vtu':
                reader = vtk.vtkXMLUnstructuredGridReader()
            elif file_ext == '.vts':
                reader = vtk.vtkXMLStructuredGridReader()
            elif file_ext == '.vtp':
                reader = vtk.vtkXMLPolyDataReader()
            else:
                raise ValueError(f"Unsupported file format: {file_ext}")
            
            reader.SetFileName(file_path)
            reader.Update()
            
            # Get the output
            output = reader.GetOutput()
            
            if output.GetNumberOfPoints() == 0:
                raise ValueError("VTK file contains no data points")
            
            # Convert to unstructured grid if needed
            if isinstance(output, vtk.vtkStructuredGrid):
                print("Converting StructuredGrid to UnstructuredGrid...")
                converter = vtk.vtkStructuredGridGeometryFilter()
                converter.SetInputData(output)
                converter.Update()
                
                # Create unstructured grid
                self.vtk_data = vtk.vtkUnstructuredGrid()
                geometry = converter.GetOutput()
                
                # Convert polydata to unstructured grid
                append_filter = vtk.vtkAppendFilter()
                append_filter.AddInputData(geometry)
                append_filter.Update()
                self.vtk_data = append_filter.GetOutput()
                
            elif isinstance(output, vtk.vtkPolyData):
                print("Converting PolyData to UnstructuredGrid...")
                append_filter = vtk.vtkAppendFilter()
                append_filter.AddInputData(output)
                append_filter.Update()
                self.vtk_data = append_filter.GetOutput()
                
            else:
                self.vtk_data = output
                
            # Verify we have scalar data
            scalar_array = self.vtk_data.GetPointData().GetScalars()
            if not scalar_array:
                # Try to get the first array
                point_data = self.vtk_data.GetPointData()
                if point_data.GetNumberOfArrays() > 0:
                    scalar_array = point_data.GetArray(0)
                    point_data.SetScalars(scalar_array)
                    print(f"Using array '{scalar_array.GetName()}' as scalar field")
                else:
                    raise ValueError("VTK file contains no scalar data")
                
            # Setup field visualization
            self.setup_field_visualization()
            
            # Update analyzer
            self.flux_analyzer.set_vtk_data(self.vtk_data)
            
            self.status_label.setText(
                f"Loaded VTK data: {self.vtk_data.GetNumberOfPoints()} points, "
                f"{self.vtk_data.GetNumberOfCells()} cells\n"
                f"Scalar field: {scalar_array.GetName() if scalar_array else 'None'}"
            )
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to load VTK data:\n{str(e)}")
            import traceback
            print(f"VTK loading error details:\n{traceback.format_exc()}")
            
    def create_legacy_vtk_reader(self, file_path):
        """Create appropriate reader for legacy VTK files"""
        # Read first few lines to determine type
        with open(file_path, 'r') as f:
            lines = [f.readline().strip() for _ in range(10)]
            
        content = ' '.join(lines).upper()
        
        if 'STRUCTURED_GRID' in content:
            print("Detected StructuredGrid format")
            return vtk.vtkStructuredGridReader()
        elif 'UNSTRUCTURED_GRID' in content:
            print("Detected UnstructuredGrid format")
            return vtk.vtkUnstructuredGridReader()
        elif 'POLYDATA' in content:
            print("Detected PolyData format")
            return vtk.vtkPolyDataReader()
        elif 'STRUCTURED_POINTS' in content:
            print("Detected StructuredPoints format")
            return vtk.vtkStructuredPointsReader()
        else:
            print("Unknown format, trying UnstructuredGridReader")
            return vtk.vtkUnstructuredGridReader()
        
    # def load_vtk_data(self):
    #     """Load VTK data file"""
    #     file_path, _ = QFileDialog.getOpenFileName(
    #         self, "Load VTK Data", "", "VTK Files (*.vtk *.vtu *.vtp);;All Files (*)"
    #     )
        
    #     if not file_path:
    #         return
            
    #     try:
    #         # Read VTK file
    #         reader = vtk.vtkUnstructuredGridReader()
    #         reader.SetFileName(file_path)
    #         reader.Update()
            
    #         self.vtk_data = reader.GetOutput()
            
    #         if self.vtk_data.GetNumberOfPoints() == 0:
    #             raise ValueError("VTK file contains no data points")
                
    #         # Setup field visualization
    #         self.setup_field_visualization()
            
    #         # Update analyzer
    #         self.flux_analyzer.set_vtk_data(self.vtk_data)
            
    #         self.status_label.setText(
    #             f"Loaded VTK data: {self.vtk_data.GetNumberOfPoints()} points, "
    #             f"{self.vtk_data.GetNumberOfCells()} cells"
    #         )
            
    #     except Exception as e:
    #         QMessageBox.warning(self, "Error", f"Failed to load VTK data:\n{str(e)}")
            
    def setup_field_visualization(self):
        """Setup 3D field visualization"""
        if not self.vtk_data:
            return
            
        # Create field mapper with transparency
        field_mapper = vtk.vtkDataSetMapper()
        field_mapper.SetInputData(self.vtk_data)
        field_mapper.SetScalarModeToUsePointData()
        
        # Setup lookup table for coloring
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.667, 0.0)  # Blue to red
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(1.0, 1.0)
        lut.SetAlphaRange(0.1, 0.8)
        lut.SetNumberOfTableValues(256)
        lut.Build()
        
        field_mapper.SetLookupTable(lut)
        field_mapper.SetScalarRange(self.vtk_data.GetScalarRange())
        
        # Create field actor
        if hasattr(self, 'field_actor'):
            self.renderer.RemoveActor(self.field_actor)
            
        self.field_actor = vtk.vtkActor()
        self.field_actor.SetMapper(field_mapper)
        self.field_actor.GetProperty().SetOpacity(0.3)  # Translucent
        
        self.renderer.AddActor(self.field_actor)
        
        # Add scalar bar
        if hasattr(self, 'scalar_bar'):
            self.renderer.RemoveActor2D(self.scalar_bar)
            
        self.scalar_bar = vtk.vtkScalarBarActor()
        self.scalar_bar.SetLookupTable(lut)
        self.scalar_bar.SetTitle("Electron Flux")
        self.scalar_bar.SetPosition(0.85, 0.1)
        self.scalar_bar.SetWidth(0.1)
        self.scalar_bar.SetHeight(0.8)
        
        self.renderer.AddActor2D(self.scalar_bar)
        
        self.vtk_widget.GetRenderWindow().Render()
        
    def load_orbital_data(self):
        """Load orbital CSV data"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Orbital Data", "", "CSV Files (*.csv);;All Files (*)"
        )
        
        if not file_path:
            return
            
        try:
            # Read CSV file
            df = pd.read_csv(file_path)
            
            # Check for required columns
            required_cols = ['time', 'x', 'y', 'z']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")
                
            # Convert to orbital points
            self.orbital_path = []
            for _, row in df.iterrows():
                vx = row.get('vx', 0)
                vy = row.get('vy', 0) 
                vz = row.get('vz', 0)
                
                point = OrbitalPoint(
                    time=row['time'], x=row['x'], y=row['y'], z=row['z'],
                    vx=vx, vy=vy, vz=vz
                )
                self.orbital_path.append(point)
                
            # Update time slider
            self.time_slider.setMaximum(len(self.orbital_path) - 1)
            
            # Update analyzer
            self.flux_analyzer.set_orbital_data(self.orbital_path)
            
            # Create path visualization
            self.create_path_visualization()
            
            # Reset animation
            self.reset_animation()
            
            self.status_label.setText(
                f"Loaded orbital data: {len(self.orbital_path)} points "
                f"over {self.orbital_path[-1].time:.2f} hours"
            )
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to load orbital data:\n{str(e)}")
            
    def create_path_visualization(self):
        """Create 3D orbital path visualization"""
        if not self.orbital_path:
            return
            
        # Create path polyline
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        
        for i, point in enumerate(self.orbital_path):
            points.InsertNextPoint(point.x, point.y, point.z)
            
        # Create line segments
        for i in range(len(self.orbital_path) - 1):
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, i)
            line.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(line)
            
        path_polydata = vtk.vtkPolyData()
        path_polydata.SetPoints(points)
        path_polydata.SetLines(lines)
        
        path_mapper = vtk.vtkPolyDataMapper()
        path_mapper.SetInputData(path_polydata)
        
        if hasattr(self, 'path_actor'):
            self.renderer.RemoveActor(self.path_actor)
            
        self.path_actor = vtk.vtkActor()
        self.path_actor.SetMapper(path_mapper)
        self.path_actor.GetProperty().SetColor(0.8, 0.8, 0.0)  # Yellow
        self.path_actor.GetProperty().SetLineWidth(2.0)
        self.path_actor.GetProperty().SetOpacity(0.7)
        
        self.renderer.AddActor(self.path_actor)
        
        # Create object representation
        self.create_object_representation()
        
        self.vtk_widget.GetRenderWindow().Render()
        
    def create_object_representation(self):
        """Create representation of the orbiting object"""
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(self.cross_section_spinbox.value() * 10.0)  # Scale for visibility
        
        object_mapper = vtk.vtkPolyDataMapper()
        object_mapper.SetInputConnection(sphere.GetOutputPort())
        
        if hasattr(self, 'object_actor'):
            self.renderer.RemoveActor(self.object_actor)
            
        self.object_actor = vtk.vtkActor()
        self.object_actor.SetMapper(object_mapper)
        self.object_actor.GetProperty().SetColor(1.0, 0.2, 0.2)  # Red
        self.object_actor.GetProperty().SetOpacity(0.9)
        
        self.renderer.AddActor(self.object_actor)
        
    def start_animation(self):
        """Start the orbital animation"""
        if not self.orbital_path or not self.vtk_data:
            QMessageBox.warning(self, "Warning", "Please load both VTK data and orbital path first.")
            return
            
        self.is_playing = True
        self.animation_timer.start(self.speed_spinbox.value())
        
        # Show plot windows
        if self.slice_window:
            self.slice_window.show()
        if self.spectrum_window:
            self.spectrum_window.show()
        if self.flux_time_window:
            self.flux_time_window.show()
            
        self.play_button.setEnabled(False)
        self.pause_button.setEnabled(True)
        self.stop_button.setEnabled(True)
        
    def pause_animation(self):
        """Pause the animation"""
        self.is_playing = False
        self.animation_timer.stop()
        
        self.play_button.setEnabled(True)
        self.pause_button.setEnabled(False)
        
    def stop_animation(self):
        """Stop the animation"""
        self.is_playing = False
        self.animation_timer.stop()
        self.reset_animation()
        
        self.play_button.setEnabled(True)
        self.pause_button.setEnabled(False)
        self.stop_button.setEnabled(False)
        
    def animation_step(self):
        """Perform one animation step"""
        if self.current_time_index >= len(self.orbital_path) - 1:
            # Completed one orbit
            self.stop_animation()
            return
            
        self.current_time_index += 1
        self.update_visualization()
        self.update_plots()
        
    def update_visualization(self):
        """Update the 3D visualization"""
        if not self.orbital_path or self.current_time_index >= len(self.orbital_path):
            return
            
        current_point = self.orbital_path[self.current_time_index]
        
        # Update object position
        if hasattr(self, 'object_actor'):
            self.object_actor.SetPosition(current_point.x, current_point.y, current_point.z)
            
        # Update time slider
        self.time_slider.setValue(self.current_time_index)
        
        # Update time label
        hours = current_point.time
        h = int(hours)
        m = int((hours - h) * 60)
        s = int(((hours - h) * 60 - m) * 60)
        self.time_label.setText(f"{h:02d}:{m:02d}:{s:02d}")
        
        # Update status
        flux = self.flux_analyzer.analyze_flux_at_point(current_point)
        self.status_label.setText(
            f"Time: {current_point.time:.2f}h | "
            f"Position: ({current_point.x:.1f}, {current_point.y:.1f}, {current_point.z:.1f}) km | "
            f"Flux: {flux:.2e} particles/s"
        )
        
        self.vtk_widget.GetRenderWindow().Render()
        
    def update_plots(self):
        """Update all plot windows"""
        if not self.orbital_path or self.current_time_index >= len(self.orbital_path):
            return
            
        current_point = self.orbital_path[self.current_time_index]
        
        # Update slice plot
        if self.slice_window and self.vtk_data:
            self.slice_window.update_slice(current_point.phi, self.vtk_data)
            self.slice_window.set_object_position(current_point.x, current_point.y, current_point.z)
            
        # Update spectrum plot (placeholder)
        if self.spectrum_window:
            # Generate mock spectrum data
            energies = np.linspace(0.1, 5.0, 50)
            spectrum = np.exp(-energies / 0.5) * (1 + 0.2 * np.sin(energies))
            self.spectrum_window.update_spectrum(energies, spectrum)
            
        # Update flux time plot
        if self.flux_time_window:
            flux = self.flux_analyzer.analyze_flux_at_point(current_point)
            self.flux_time_window.add_flux_point(current_point.time, flux)
            
    def jump_to_time(self, time_index):
        """Jump to specific time index"""
        if 0 <= time_index < len(self.orbital_path):
            self.current_time_index = time_index
            self.update_visualization()
            self.update_plots()
            
    def set_animation_speed(self, speed_ms):
        """Set animation speed"""
        if self.is_playing:
            self.animation_timer.setInterval(speed_ms)
            
    def update_cross_section(self, radius):
        """Update object cross section"""
        self.flux_analyzer.set_cross_section(radius)
        if hasattr(self, 'object_actor'):
            self.renderer.RemoveActor(self.object_actor)
            self.create_object_representation()
            self.vtk_widget.GetRenderWindow().Render()
            
    def reset_animation(self):
        """Reset animation to beginning"""
        self.current_time_index = 0
        if self.flux_time_window:
            self.flux_time_window.clear_data()
        self.update_visualization()
        
    def show_slice_window(self):
        """Show the slice plot window"""
        if not self.slice_window:
            self.slice_window = SlicePlotWindow(self)
        self.slice_window.show()
        self.slice_window.raise_()
        
    def show_spectrum_window(self):
        """Show the spectrum plot window"""
        if not self.spectrum_window:
            self.spectrum_window = SpectrumPlotWindow(self)
        self.spectrum_window.show()
        self.spectrum_window.raise_()
        
    def show_flux_time_window(self):
        """Show the flux vs time window"""
        if not self.flux_time_window:
            self.flux_time_window = FluxTimePlotWindow(self)
        self.flux_time_window.show()
        self.flux_time_window.raise_()

def main():
    """Main application entry point"""
    app = QApplication(sys.argv)
    
    # Set application properties
    app.setApplicationName("Electron Flux Visualizer")
    app.setApplicationVersion("1.0")
    app.setOrganizationName("Scientific Visualization Lab")
    
    # Create and show main window
    window = ElectronFluxVisualizer()
    window.show()
    
    return app.exec()

if __name__ == '__main__':
    sys.exit(main())
