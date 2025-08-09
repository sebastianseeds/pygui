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
    QSplitter, QFrame, QComboBox, QCheckBox
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

        self.setup_visualization_controls()
        
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
        self.control_layout = QVBoxLayout(control_panel)
        
        # File loading group
        file_group = QGroupBox("Data Loading")
        file_layout = QFormLayout(file_group)
        
        self.load_vtk_button = QPushButton("Load VTK Data")
        self.load_orbit_button = QPushButton("Load Orbital CSV")
        
        file_layout.addRow(self.load_vtk_button)
        file_layout.addRow(self.load_orbit_button)
        
        self.control_layout.addWidget(file_group)
        
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
        
        self.control_layout.addWidget(anim_group)
        
        # Analysis parameters group
        analysis_group = QGroupBox("Analysis Parameters")
        analysis_layout = QFormLayout(analysis_group)
        
        self.cross_section_spinbox = QDoubleSpinBox()
        self.cross_section_spinbox.setRange(0.1, 100.0)
        self.cross_section_spinbox.setValue(1.0)
        self.cross_section_spinbox.setSingleStep(0.1)
        self.cross_section_spinbox.setSuffix(" m")
        
        analysis_layout.addRow("Cross Section Radius:", self.cross_section_spinbox)
        self.control_layout.addWidget(analysis_group)
        
        # Plot windows group
        plots_group = QGroupBox("Plot Windows")
        plots_layout = QVBoxLayout(plots_group)
        
        self.show_slice_button = QPushButton("Show Slice Plot")
        self.show_spectrum_button = QPushButton("Show Spectrum Plot")
        self.show_flux_time_button = QPushButton("Show Flux vs Time")
        
        plots_layout.addWidget(self.show_slice_button)
        plots_layout.addWidget(self.show_spectrum_button)
        plots_layout.addWidget(self.show_flux_time_button)
        
        self.control_layout.addWidget(plots_group)
        
        # Status
        self.status_label = QLabel("Ready - Load VTK data and orbital CSV to begin")
        self.status_label.setWordWrap(True)
        self.control_layout.addWidget(self.status_label)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.control_layout.addWidget(self.progress_bar)
        
        self.control_layout.addStretch()
        
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

    def setup_visualization_controls(self):
        """Add ParaView-inspired visualization controls to the UI"""
        
        # Add to the control panel after analysis parameters
        viz_group = QGroupBox("Field Visualization")
        viz_layout = QVBoxLayout(viz_group)
        
        # Visualization mode selection
        mode_layout = QFormLayout()
        self.viz_mode_combo = QComboBox()
        self.viz_mode_combo.addItems([
            "Volume Rendering",
            "Isosurfaces", 
            "Point Cloud",
            "Wireframe",
            "Surface with Edges",
            "Slice Planes"
        ])
        self.viz_mode_combo.currentTextChanged.connect(self.change_visualization_mode)
        mode_layout.addRow("Visualization Mode:", self.viz_mode_combo)
        
        # Opacity control
        self.opacity_slider = QSlider(Qt.Orientation.Horizontal)
        self.opacity_slider.setRange(0, 100)
        self.opacity_slider.setValue(40)
        self.opacity_slider.valueChanged.connect(self.update_opacity)
        self.opacity_label = QLabel("40%")
        opacity_layout = QHBoxLayout()
        opacity_layout.addWidget(self.opacity_slider)
        opacity_layout.addWidget(self.opacity_label)
        mode_layout.addRow("Opacity:", opacity_layout)
        
        # Threshold controls
        self.threshold_enabled = QCheckBox("Enable Thresholding")
        self.threshold_enabled.toggled.connect(self.toggle_threshold)
        
        self.threshold_min_slider = QSlider(Qt.Orientation.Horizontal)
        self.threshold_min_slider.setRange(0, 100)
        self.threshold_min_slider.setValue(10)
        self.threshold_min_slider.valueChanged.connect(self.update_threshold)
        
        self.threshold_max_slider = QSlider(Qt.Orientation.Horizontal) 
        self.threshold_max_slider.setRange(0, 100)
        self.threshold_max_slider.setValue(90)
        self.threshold_max_slider.valueChanged.connect(self.update_threshold)
        
        threshold_layout = QFormLayout()
        threshold_layout.addRow(self.threshold_enabled)
        threshold_layout.addRow("Min Threshold:", self.threshold_min_slider)
        threshold_layout.addRow("Max Threshold:", self.threshold_max_slider)
        
        # Isosurface controls
        self.isosurface_enabled = QCheckBox("Show Isosurfaces")
        self.isosurface_enabled.toggled.connect(self.toggle_isosurfaces)
        
        self.num_isosurfaces_spin = QSpinBox()
        self.num_isosurfaces_spin.setRange(1, 10)
        self.num_isosurfaces_spin.setValue(3)
        self.num_isosurfaces_spin.valueChanged.connect(self.update_isosurfaces)
        
        iso_layout = QFormLayout()
        iso_layout.addRow(self.isosurface_enabled)
        iso_layout.addRow("Number of Isosurfaces:", self.num_isosurfaces_spin)
        
        # Color map selection
        self.colormap_combo = QComboBox()
        self.colormap_combo.addItems([
            "Blue to Red",
            "Viridis", 
            "Plasma",
            "Cool to Warm",
            "Rainbow",
            "Grayscale"
        ])
        self.colormap_combo.currentTextChanged.connect(self.change_colormap)
        mode_layout.addRow("Color Map:", self.colormap_combo)
        
        # Add all layouts to group
        viz_layout.addLayout(mode_layout)
        viz_layout.addLayout(threshold_layout)
        viz_layout.addLayout(iso_layout)
        
        # Insert into main control layout (after analysis_group)
        control_layout = self.findChild(QVBoxLayout)  # Find the main control layout
        # Insert before the plots group
        plots_group_index = None
        for i in range(control_layout.count()):
            widget = control_layout.itemAt(i).widget()
            if isinstance(widget, QGroupBox) and widget.title() == "Plot Windows":
                plots_group_index = i
                break
        
        if plots_group_index:
            control_layout.insertWidget(plots_group_index, viz_group)
        else:
            control_layout.addWidget(viz_group)

    def change_visualization_mode(self, mode):
        """Change the field visualization mode - IMPROVED ERROR HANDLING"""
        if not self.vtk_data:
            print("No VTK data loaded")
            return
            
        print(f"Changing visualization mode to: {mode}")
        
        try:
            # Remove existing actors/volumes
            if hasattr(self, 'field_actor') and self.field_actor:
                self.renderer.RemoveActor(self.field_actor)
                self.field_actor = None
                
            if hasattr(self, 'volume_actor') and self.volume_actor:
                self.renderer.RemoveVolume(self.volume_actor)
                self.volume_actor = None
                
            if hasattr(self, 'slice_actors'):
                for actor in self.slice_actors:
                    self.renderer.RemoveActor(actor)
                self.slice_actors = []
                
            # Apply new visualization mode
            if mode == "Volume Rendering":
                self.setup_volume_rendering()
            elif mode == "Isosurfaces":
                self.setup_isosurface_rendering()
            elif mode == "Point Cloud":
                self.setup_point_cloud_rendering()
            elif mode == "Wireframe":
                self.setup_wireframe_rendering()
            elif mode == "Surface with Edges":
                self.setup_surface_with_edges()
            elif mode == "Slice Planes":
                self.setup_slice_planes()
            else:
                # Fallback to original visualization
                self.setup_field_visualization()
                
            self.vtk_widget.GetRenderWindow().Render()
            print(f"Visualization mode changed to: {mode}")
            
        except Exception as e:
            print(f"Error changing visualization mode: {e}")
            # Fallback to basic field visualization
            self.setup_field_visualization()
            self.vtk_widget.GetRenderWindow().Render()

    def setup_volume_rendering(self):
        """Setup volume rendering for 3D scalar field - FIXED VERSION"""
        if not self.vtk_data:
            return
            
        print("Setting up volume rendering...")
        
        try:
            # For volume rendering, we need image data
            if not isinstance(self.vtk_data, vtk.vtkImageData):
                # Create a volume from point cloud using vtkGaussianSplatter
                splatter = vtk.vtkGaussianSplatter()
                splatter.SetInputData(self.vtk_data)
                splatter.SetSampleDimensions(40, 40, 40)  # Reduced for performance
                
                # Calculate bounds for radius
                bounds = self.vtk_data.GetBounds()
                max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
                radius = max_extent / 20.0  # Adaptive radius
                
                splatter.SetRadius(radius)
                splatter.SetExponentFactor(-5)  # Controls falloff
                splatter.Update()
                volume_data = splatter.GetOutput()
                
                if volume_data.GetNumberOfPoints() == 0:
                    print("Warning: Volume splatter produced no data")
                    return
            else:
                volume_data = self.vtk_data
                
            # Volume mapper - try GPU first, fallback to CPU
            try:
                volume_mapper = vtk.vtkGPUVolumeRayCastMapper()
                print("Using GPU volume mapper")
            except:
                volume_mapper = vtk.vtkFixedPointVolumeRayCastMapper()
                print("Using CPU volume mapper")
                
            volume_mapper.SetInputData(volume_data)
            
            # Volume properties
            volume_property = vtk.vtkVolumeProperty()
            volume_property.ShadeOn()
            volume_property.SetInterpolationTypeToLinear()
            
            # Color transfer function
            color_func = vtk.vtkColorTransferFunction()
            scalar_range = volume_data.GetScalarRange()
            
            if scalar_range[1] > scalar_range[0]:  # Valid range
                color_func.AddRGBPoint(scalar_range[0], 0.0, 0.0, 0.2)  # Dark blue
                color_func.AddRGBPoint(scalar_range[1] * 0.3, 0.0, 0.0, 1.0)  # Blue
                color_func.AddRGBPoint(scalar_range[1] * 0.6, 0.0, 1.0, 0.0)  # Green
                color_func.AddRGBPoint(scalar_range[1] * 0.9, 1.0, 1.0, 0.0)  # Yellow
                color_func.AddRGBPoint(scalar_range[1], 1.0, 0.0, 0.0)  # Red
            
            # Opacity transfer function
            opacity_func = vtk.vtkPiecewiseFunction()
            if scalar_range[1] > scalar_range[0]:
                opacity_func.AddPoint(scalar_range[0], 0.0)
                opacity_func.AddPoint(scalar_range[1] * 0.1, 0.0)
                opacity_func.AddPoint(scalar_range[1] * 0.3, 0.1)
                opacity_func.AddPoint(scalar_range[1] * 0.6, 0.3)
                opacity_func.AddPoint(scalar_range[1], 0.8)
            
            volume_property.SetColor(color_func)
            volume_property.SetScalarOpacity(opacity_func)
            
            # Create volume
            if hasattr(self, 'volume_actor'):
                self.renderer.RemoveVolume(self.volume_actor)
                
            self.volume_actor = vtk.vtkVolume()
            self.volume_actor.SetMapper(volume_mapper)
            self.volume_actor.SetProperty(volume_property)
            
            self.renderer.AddVolume(self.volume_actor)
            
            print("Volume rendering setup complete")
            
        except Exception as e:
            print(f"Volume rendering failed: {e}")
            # Fallback to point cloud
            self.setup_point_cloud_rendering()

    def setup_isosurface_rendering(self):
        """Setup isosurface rendering with multiple contour levels - FIXED VERSION"""
        if not self.vtk_data:
            return
            
        print("Setting up isosurface rendering...")
        
        try:
            # Create contour filter
            contour = vtk.vtkContourFilter()
            contour.SetInputData(self.vtk_data)
            
            # Generate isosurface values
            scalar_range = self.vtk_data.GetScalarRange()
            num_contours = self.num_isosurfaces_spin.value()
            
            if scalar_range[1] > scalar_range[0]:
                contour.GenerateValues(num_contours, 
                                      scalar_range[1] * 0.2,  # Start at 20% of max
                                      scalar_range[1] * 0.9)  # End at 90% of max
            contour.Update()
            
            if contour.GetOutput().GetNumberOfPoints() == 0:
                print("Warning: No isosurfaces generated")
                return
            
            # Create mapper and actor
            contour_mapper = vtk.vtkPolyDataMapper()
            contour_mapper.SetInputConnection(contour.GetOutputPort())
            contour_mapper.SetScalarRange(scalar_range)
            
            # Setup lookup table
            lut = self.create_lookup_table(self.colormap_combo.currentText())
            contour_mapper.SetLookupTable(lut)
            
            if hasattr(self, 'field_actor'):
                self.renderer.RemoveActor(self.field_actor)
                
            self.field_actor = vtk.vtkActor()
            self.field_actor.SetMapper(contour_mapper)
            self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
            
            self.renderer.AddActor(self.field_actor)
            
            print(f"Created {num_contours} isosurfaces")
            
        except Exception as e:
            print(f"Isosurface rendering failed: {e}")
            # Fallback to point cloud
            self.setup_point_cloud_rendering()

    def setup_point_cloud_rendering(self):
        """Setup point cloud rendering with spherical glyphs"""
        if not self.vtk_data:
            return
            
        print("Setting up point cloud rendering...")
        
        # Create sphere glyphs
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(100)  # 100 km radius
        sphere.SetThetaResolution(8)
        sphere.SetPhiResolution(8)
        
        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(self.vtk_data)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.SetScaleModeToScaleByScalar()
        glyph.SetScaleFactor(0.1)
        glyph.Update()
        
        # Create mapper
        glyph_mapper = vtk.vtkPolyDataMapper()
        glyph_mapper.SetInputConnection(glyph.GetOutputPort())
        glyph_mapper.SetScalarRange(self.vtk_data.GetScalarRange())
        
        # Setup lookup table
        lut = self.create_lookup_table(self.colormap_combo.currentText())
        glyph_mapper.SetLookupTable(lut)
        
        self.field_actor = vtk.vtkActor()
        self.field_actor.SetMapper(glyph_mapper)
        self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
        
        self.renderer.AddActor(self.field_actor)

    def setup_wireframe_rendering(self):
        """Setup wireframe rendering"""
        if not self.vtk_data:
            return
            
        print("Setting up wireframe rendering...")
        
        # Extract surface
        surface = vtk.vtkDataSetSurfaceFilter()
        surface.SetInputData(self.vtk_data)
        surface.Update()
        
        # Create mapper
        wireframe_mapper = vtk.vtkPolyDataMapper()
        wireframe_mapper.SetInputConnection(surface.GetOutputPort())
        wireframe_mapper.SetScalarRange(self.vtk_data.GetScalarRange())
        
        # Setup lookup table
        lut = self.create_lookup_table(self.colormap_combo.currentText())
        wireframe_mapper.SetLookupTable(lut)
        
        self.field_actor = vtk.vtkActor()
        self.field_actor.SetMapper(wireframe_mapper)
        self.field_actor.GetProperty().SetRepresentationToWireframe()
        self.field_actor.GetProperty().SetLineWidth(1.0)
        self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
        
        self.renderer.AddActor(self.field_actor)

    def setup_surface_with_edges(self):
        """Setup surface rendering with edge display"""
        if not self.vtk_data:
            return
            
        print("Setting up surface with edges...")
        
        # Extract surface
        surface = vtk.vtkDataSetSurfaceFilter()
        surface.SetInputData(self.vtk_data)
        surface.Update()
        
        # Create mapper
        surface_mapper = vtk.vtkPolyDataMapper()
        surface_mapper.SetInputConnection(surface.GetOutputPort())
        surface_mapper.SetScalarRange(self.vtk_data.GetScalarRange())
        
        # Setup lookup table
        lut = self.create_lookup_table(self.colormap_combo.currentText())
        surface_mapper.SetLookupTable(lut)
        
        self.field_actor = vtk.vtkActor()
        self.field_actor.SetMapper(surface_mapper)
        self.field_actor.GetProperty().SetRepresentationToSurface()
        self.field_actor.GetProperty().EdgeVisibilityOn()
        self.field_actor.GetProperty().SetEdgeColor(0.2, 0.2, 0.2)
        self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
        
        self.renderer.AddActor(self.field_actor)

    def setup_slice_planes(self):
        """Setup slice plane visualization"""
        if not self.vtk_data:
            return
            
        print("Setting up slice planes...")
        
        # Create three orthogonal slice planes
        planes = []
        
        # XY plane (Z = 0)
        plane_xy = vtk.vtkPlane()
        plane_xy.SetOrigin(0, 0, 0)
        plane_xy.SetNormal(0, 0, 1)
        planes.append(plane_xy)
        
        # XZ plane (Y = 0) 
        plane_xz = vtk.vtkPlane()
        plane_xz.SetOrigin(0, 0, 0)
        plane_xz.SetNormal(0, 1, 0)
        planes.append(plane_xz)
        
        # YZ plane (X = 0)
        plane_yz = vtk.vtkPlane()
        plane_yz.SetOrigin(0, 0, 0)
        plane_yz.SetNormal(1, 0, 0)
        planes.append(plane_yz)
        
        self.slice_actors = []
        
        for i, plane in enumerate(planes):
            # Create cutter
            cutter = vtk.vtkCutter()
            cutter.SetInputData(self.vtk_data)
            cutter.SetCutFunction(plane)
            cutter.Update()
            
            # Create mapper
            slice_mapper = vtk.vtkPolyDataMapper()
            slice_mapper.SetInputConnection(cutter.GetOutputPort())
            slice_mapper.SetScalarRange(self.vtk_data.GetScalarRange())
            
            # Setup lookup table
            lut = self.create_lookup_table(self.colormap_combo.currentText())
            slice_mapper.SetLookupTable(lut)
            
            # Create actor
            slice_actor = vtk.vtkActor()
            slice_actor.SetMapper(slice_mapper)
            slice_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
            
            self.slice_actors.append(slice_actor)
            self.renderer.AddActor(slice_actor)

    def create_lookup_table(self, colormap_name):
        """Create lookup table based on colormap name"""
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(256)
        
        if colormap_name == "Blue to Red":
            lut.SetHueRange(0.667, 0.0)  # Blue to red
        elif colormap_name == "Viridis":
            # Approximate viridis colormap
            lut.SetHueRange(0.7, 0.15)
            lut.SetSaturationRange(1.0, 0.7)
        elif colormap_name == "Plasma":
            # Approximate plasma colormap  
            lut.SetHueRange(0.8, 0.05)
            lut.SetSaturationRange(1.0, 1.0)
        elif colormap_name == "Cool to Warm":
            lut.SetHueRange(0.667, 0.0)
            lut.SetSaturationRange(1.0, 1.0)
        elif colormap_name == "Rainbow":
            lut.SetHueRange(0.0, 0.667)
        elif colormap_name == "Grayscale":
            lut.SetHueRange(0.0, 0.0)
            lut.SetSaturationRange(0.0, 0.0)
            
        lut.SetValueRange(0.0, 1.0)
        lut.Build()
        return lut

    def update_opacity(self, value):
        """Update field visualization opacity"""
        opacity = value / 100.0
        self.opacity_label.setText(f"{value}%")
        
        if hasattr(self, 'field_actor') and self.field_actor:
            self.field_actor.GetProperty().SetOpacity(opacity)
            
        if hasattr(self, 'volume_actor') and self.volume_actor:
            # Update volume opacity (more complex)
            volume_property = self.volume_actor.GetProperty()
            opacity_func = volume_property.GetScalarOpacity()
            # Scale existing opacity function
            
        if hasattr(self, 'slice_actors'):
            for actor in self.slice_actors:
                actor.GetProperty().SetOpacity(opacity)
                
        self.vtk_widget.GetRenderWindow().Render()

    def toggle_threshold(self, enabled):
        """Toggle threshold filtering"""
        if enabled:
            self.update_threshold()
        else:
            # Remove threshold and restore full data
            self.change_visualization_mode(self.viz_mode_combo.currentText())

    def update_threshold(self):
        """Update threshold values with correct VTK API"""
        if not self.threshold_enabled.isChecked() or not self.vtk_data:
            return
            
        scalar_range = self.vtk_data.GetScalarRange()
        min_val = scalar_range[0] + (scalar_range[1] - scalar_range[0]) * (self.threshold_min_slider.value() / 100.0)
        max_val = scalar_range[0] + (scalar_range[1] - scalar_range[0]) * (self.threshold_max_slider.value() / 100.0)
        
        print(f"Applying threshold: {min_val:.2e} to {max_val:.2e}")
        
        # Apply threshold filter with correct VTK 9+ API
        threshold = vtk.vtkThreshold()
        threshold.SetInputData(self.vtk_data)
        
        # Use the correct modern VTK threshold API
        threshold.SetLowerThreshold(min_val)
        threshold.SetUpperThreshold(max_val)
        threshold.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)
        threshold.Update()
        
        # Store the thresholded data temporarily
        self.thresholded_data = threshold.GetOutput()
        
        # Re-apply current visualization mode with thresholded data
        self.apply_visualization_to_data(self.thresholded_data)

    def apply_visualization_to_data(self, data):
        """Apply current visualization mode to given data"""
        if not data:
            return
            
        # Remove existing actors
        if hasattr(self, 'field_actor') and self.field_actor:
            self.renderer.RemoveActor(self.field_actor)
        if hasattr(self, 'volume_actor') and self.volume_actor:
            self.renderer.RemoveVolume(self.volume_actor)
        if hasattr(self, 'slice_actors'):
            for actor in self.slice_actors:
                self.renderer.RemoveActor(actor)
            
        # Temporarily replace vtk_data
        original_data = self.vtk_data
        self.vtk_data = data
        
        # Apply current visualization mode
        current_mode = self.viz_mode_combo.currentText()
        
        if current_mode == "Volume Rendering":
            self.setup_volume_rendering()
        elif current_mode == "Isosurfaces":
            self.setup_isosurface_rendering()
        elif current_mode == "Point Cloud":
            self.setup_point_cloud_rendering()
        elif current_mode == "Wireframe":
            self.setup_wireframe_rendering()
        elif current_mode == "Surface with Edges":
            self.setup_surface_with_edges()
        elif current_mode == "Slice Planes":
            self.setup_slice_planes()
        else:
            # Default to original field visualization
            self.setup_field_visualization()
            
        # Restore original data reference
        self.vtk_data = original_data
        
        self.vtk_widget.GetRenderWindow().Render()

    def toggle_isosurfaces(self, enabled):
        """Toggle isosurface display"""
        if enabled and self.viz_mode_combo.currentText() != "Isosurfaces":
            self.viz_mode_combo.setCurrentText("Isosurfaces")
        elif not enabled and hasattr(self, 'isosurface_actors'):
            for actor in self.isosurface_actors:
                self.renderer.RemoveActor(actor)
            self.vtk_widget.GetRenderWindow().Render()

    def update_isosurfaces(self):
        """Update number of isosurfaces"""
        if (self.viz_mode_combo.currentText() == "Isosurfaces" and 
            self.isosurface_enabled.isChecked()):
            self.setup_isosurface_rendering()

    def change_colormap(self, colormap_name):
        """Change the color mapping"""
        if hasattr(self, 'field_actor') and self.field_actor:
            mapper = self.field_actor.GetMapper()
            if mapper:
                lut = self.create_lookup_table(colormap_name)
                mapper.SetLookupTable(lut)
                
                # Update scalar bar
                if hasattr(self, 'scalar_bar'):
                    self.scalar_bar.SetLookupTable(lut)
                    
        self.vtk_widget.GetRenderWindow().Render()
        
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
        """Load VTK data file with comprehensive format support"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load VTK Data", "", 
            "VTK Files (*.vtk *.vtu *.vtp *.vts *.vti);;XML VTK (*.vtu *.vtp *.vts *.vti);;Legacy VTK (*.vtk);;All Files (*)"
        )
        
        if not file_path:
            return
            
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        
        try:
            # Determine file type and create appropriate reader
            file_ext = Path(file_path).suffix.lower()
            
            print(f"Loading VTK file: {file_path}")
            print(f"File extension: {file_ext}")
            
            if file_ext == '.vts':
                print("Using XML Structured Grid Reader")
                reader = vtk.vtkXMLStructuredGridReader()
            elif file_ext == '.vtu':
                print("Using XML Unstructured Grid Reader")
                reader = vtk.vtkXMLUnstructuredGridReader()
            elif file_ext == '.vtp':
                print("Using XML PolyData Reader")
                reader = vtk.vtkXMLPolyDataReader()
            elif file_ext == '.vti':
                print("Using XML Image Data Reader")
                reader = vtk.vtkXMLImageDataReader()
            elif file_ext == '.vtk':
                # Legacy format - auto-detect type
                reader = self.create_legacy_vtk_reader(file_path)
            else:
                raise ValueError(f"Unsupported file format: {file_ext}")
            
            self.progress_bar.setValue(25)
            
            # Set filename and read
            reader.SetFileName(file_path)
            print("Reading VTK file...")
            reader.Update()
            
            self.progress_bar.setValue(50)
            
            # Get the output
            output = reader.GetOutput()
            print(f"VTK object type: {type(output).__name__}")
            print(f"Number of points: {output.GetNumberOfPoints()}")
            print(f"Number of cells: {output.GetNumberOfCells()}")
            
            if output.GetNumberOfPoints() == 0:
                raise ValueError("VTK file contains no data points")
            
            self.progress_bar.setValue(75)
            
            # Handle different VTK data types
            if isinstance(output, (vtk.vtkImageData, vtk.vtkStructuredGrid, vtk.vtkRectilinearGrid)):
                print("Converting structured data to unstructured grid...")
                self.vtk_data = self.convert_to_unstructured_grid(output)
            elif isinstance(output, vtk.vtkPolyData):
                print("Converting polydata to unstructured grid...")
                self.vtk_data = self.convert_polydata_to_unstructured_grid(output)
            elif isinstance(output, vtk.vtkUnstructuredGrid):
                print("Using unstructured grid directly")
                self.vtk_data = output
            else:
                print(f"Attempting to use {type(output).__name__} directly")
                self.vtk_data = output
                
            # Verify and setup scalar data
            self.setup_scalar_data()
            
            # Setup field visualization
            self.setup_field_visualization()
            
            # Update analyzer
            self.flux_analyzer.set_vtk_data(self.vtk_data)
            
            self.progress_bar.setValue(100)
            
            # Get scalar info for status
            scalar_array = self.vtk_data.GetPointData().GetScalars()
            scalar_name = scalar_array.GetName() if scalar_array else "None"
            scalar_range = scalar_array.GetRange() if scalar_array else (0, 0)
            
            self.status_label.setText(
                f"✓ Loaded VTK data successfully!\n"
                f"Points: {self.vtk_data.GetNumberOfPoints():,} | "
                f"Cells: {self.vtk_data.GetNumberOfCells():,}\n"
                f"Scalar: {scalar_name} | "
                f"Range: {scalar_range[0]:.2e} to {scalar_range[1]:.2e}"
            )
            
            print("VTK file loaded successfully!")
            
        except Exception as e:
            self.status_label.setText(f"Error loading VTK file: {str(e)}")
            QMessageBox.critical(self, "VTK Loading Error", 
                               f"Failed to load VTK data:\n\n{str(e)}\n\n"
                               f"Supported formats:\n"
                               f"• .vts (XML Structured Grid)\n"
                               f"• .vtu (XML Unstructured Grid)\n"
                               f"• .vtp (XML PolyData)\n"
                               f"• .vti (XML Image Data)\n"
                               f"• .vtk (Legacy VTK)")
            import traceback
            print(f"Error details:\n{traceback.format_exc()}")
            
        finally:
            self.progress_bar.setVisible(False)

    def create_legacy_vtk_reader(self, file_path):
        """Create appropriate reader for legacy VTK files with better detection"""
        try:
            # Read header to determine type
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                lines = []
                for i in range(15):  # Read more lines for better detection
                    line = f.readline()
                    if not line:
                        break
                    lines.append(line.strip().upper())
            
            content = ' '.join(lines)
            print(f"VTK file header content: {content[:200]}...")
            
            # Detect format
            if 'STRUCTURED_GRID' in content:
                print("Detected: Structured Grid")
                return vtk.vtkStructuredGridReader()
            elif 'UNSTRUCTURED_GRID' in content:
                print("Detected: Unstructured Grid")
                return vtk.vtkUnstructuredGridReader()
            elif 'POLYDATA' in content:
                print("Detected: PolyData")
                return vtk.vtkPolyDataReader()
            elif 'STRUCTURED_POINTS' in content or 'DATASET STRUCTURED_POINTS' in content:
                print("Detected: Structured Points")
                return vtk.vtkStructuredPointsReader()
            elif 'RECTILINEAR_GRID' in content:
                print("Detected: Rectilinear Grid")
                return vtk.vtkRectilinearGridReader()
            else:
                print("Unknown format, trying generic data reader")
                # Try the generic reader first
                return vtk.vtkDataSetReader()
                
        except Exception as e:
            print(f"Error detecting VTK format: {e}")
            print("Defaulting to UnstructuredGridReader")
            return vtk.vtkUnstructuredGridReader()

    def convert_to_unstructured_grid(self, structured_data):
        """Convert structured data (ImageData, StructuredGrid) to UnstructuredGrid"""
        print("Converting structured data to unstructured grid...")
        
        if isinstance(structured_data, vtk.vtkImageData):
            # Convert ImageData to UnstructuredGrid
            converter = vtk.vtkImageDataGeometryFilter()
            converter.SetInputData(structured_data)
            converter.Update()
            
            # Convert resulting polydata to unstructured grid
            return self.convert_polydata_to_unstructured_grid(converter.GetOutput())
            
        elif isinstance(structured_data, vtk.vtkStructuredGrid):
            # Convert StructuredGrid to UnstructuredGrid
            converter = vtk.vtkStructuredGridGeometryFilter()
            converter.SetInputData(structured_data)
            converter.Update()
            
            # Convert resulting polydata to unstructured grid
            return self.convert_polydata_to_unstructured_grid(converter.GetOutput())
            
        elif isinstance(structured_data, vtk.vtkRectilinearGrid):
            # Convert RectilinearGrid to UnstructuredGrid
            converter = vtk.vtkRectilinearGridGeometryFilter()
            converter.SetInputData(structured_data)
            converter.Update()
            
            return self.convert_polydata_to_unstructured_grid(converter.GetOutput())
        
        else:
            # If it's already unstructured, return as-is
            return structured_data

    def convert_polydata_to_unstructured_grid(self, polydata):
        """Convert PolyData to UnstructuredGrid"""
        print("Converting PolyData to UnstructuredGrid...")
        
        ugrid = vtk.vtkUnstructuredGrid()
        
        # Copy points
        ugrid.SetPoints(polydata.GetPoints())
        
        # Copy cells
        for i in range(polydata.GetNumberOfCells()):
            cell = polydata.GetCell(i)
            ugrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        
        # Copy point data
        ugrid.GetPointData().ShallowCopy(polydata.GetPointData())
        
        # Copy cell data
        ugrid.GetCellData().ShallowCopy(polydata.GetCellData())
        
        return ugrid

    def setup_scalar_data(self):
        """Setup and verify scalar data for visualization"""
        if not self.vtk_data:
            return
            
        point_data = self.vtk_data.GetPointData()
        
        # Check if we have scalar data
        scalar_array = point_data.GetScalars()
        
        if not scalar_array:
            print("No scalar data found, checking available arrays...")
            
            # List all available arrays
            n_arrays = point_data.GetNumberOfArrays()
            print(f"Found {n_arrays} point data arrays:")
            
            for i in range(n_arrays):
                array = point_data.GetArray(i)
                print(f"  {i}: {array.GetName()} ({array.GetNumberOfTuples()} tuples, {array.GetNumberOfComponents()} components)")
            
            if n_arrays > 0:
                # Use the first array as scalars
                first_array = point_data.GetArray(0)
                point_data.SetScalars(first_array)
                scalar_array = first_array
                print(f"Using '{first_array.GetName()}' as scalar field")
            else:
                # Create a default scalar field based on distance from origin
                print("Creating default scalar field based on distance from Earth center...")
                self.create_default_scalar_field()
                scalar_array = self.vtk_data.GetPointData().GetScalars()
        
        # Ensure scalar array has a proper name
        if scalar_array and not scalar_array.GetName():
            scalar_array.SetName("flux_field")
            
        print(f"Scalar field setup complete: {scalar_array.GetName() if scalar_array else 'None'}")

    def create_default_scalar_field(self):
        """Create a default scalar field when none exists"""
        if not self.vtk_data:
            return
            
        n_points = self.vtk_data.GetNumberOfPoints()
        scalar_values = []
        
        earth_radius = 6371.0  # km
        
        for i in range(n_points):
            point = self.vtk_data.GetPoint(i)
            x, y, z = point
            
            # Distance from Earth center
            r = np.sqrt(x*x + y*y + z*z)
            
            # Create a simple Van Allen belt-like field
            flux = 0.0
            if r > earth_radius:
                if 1.2 * earth_radius <= r <= 6 * earth_radius:
                    flux = 1e6 * np.exp(-((r - 3*earth_radius)**2) / (earth_radius)**2)
            
            scalar_values.append(flux)
        
        # Add to VTK data
        from vtk.util.numpy_support import numpy_to_vtk
        scalar_array = numpy_to_vtk(np.array(scalar_values), deep=True)
        scalar_array.SetName("generated_flux")
        self.vtk_data.GetPointData().SetScalars(scalar_array)
        
        print(f"Created default scalar field with {len(scalar_values)} values")

    def setup_field_visualization(self):
        """Setup 3D field visualization with improved handling"""
        if not self.vtk_data:
            return
            
        print("Setting up field visualization...")
        
        # Remove existing field actor
        if hasattr(self, 'field_actor'):
            self.renderer.RemoveActor(self.field_actor)
            
        # Get scalar range for proper coloring
        scalar_array = self.vtk_data.GetPointData().GetScalars()
        if not scalar_array:
            print("Warning: No scalar data for visualization")
            return
            
        scalar_range = scalar_array.GetRange()
        print(f"Scalar range: {scalar_range[0]:.2e} to {scalar_range[1]:.2e}")
        
        # Create field mapper
        field_mapper = vtk.vtkDataSetMapper()
        field_mapper.SetInputData(self.vtk_data)
        field_mapper.SetScalarModeToUsePointData()
        
        # Setup lookup table for nice coloring
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.667, 0.0)  # Blue to red
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(0.8, 1.0)
        lut.SetAlphaRange(0.2, 0.8)  # Transparency
        lut.SetNumberOfTableValues(256)
        lut.Build()
        
        field_mapper.SetLookupTable(lut)
        field_mapper.SetScalarRange(scalar_range)
        
        # Create field actor
        self.field_actor = vtk.vtkActor()
        self.field_actor.SetMapper(field_mapper)
        self.field_actor.GetProperty().SetOpacity(0.4)  # Semi-transparent
        
        # For point clouds, make points larger
        if self.vtk_data.GetNumberOfCells() == self.vtk_data.GetNumberOfPoints():
            # Likely a point cloud
            self.field_actor.GetProperty().SetPointSize(3.0)
            self.field_actor.GetProperty().SetRenderPointsAsSpheres(True)
        
        self.renderer.AddActor(self.field_actor)
        
        # Setup scalar bar
        if hasattr(self, 'scalar_bar'):
            self.renderer.RemoveActor2D(self.scalar_bar)
            
        self.scalar_bar = vtk.vtkScalarBarActor()
        self.scalar_bar.SetLookupTable(lut)
        self.scalar_bar.SetTitle(scalar_array.GetName())
        self.scalar_bar.SetPosition(0.85, 0.1)
        self.scalar_bar.SetWidth(0.12)
        self.scalar_bar.SetHeight(0.8)
        
        # Improve scalar bar appearance
        self.scalar_bar.SetNumberOfLabels(8)
        self.scalar_bar.GetLabelTextProperty().SetColor(1, 1, 1)
        self.scalar_bar.GetTitleTextProperty().SetColor(1, 1, 1)
        self.scalar_bar.GetTitleTextProperty().SetFontSize(12)
        
        self.renderer.AddActor2D(self.scalar_bar)
        
        # Reset camera to show data
        self.renderer.ResetCamera()
        
        # Render
        self.vtk_widget.GetRenderWindow().Render()
        
        print("Field visualization setup complete")
        
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
        """Create a more visible representation of the orbiting object"""
        # Create a larger, more visible sphere
        sphere = vtk.vtkSphereSource()
        # Make it much larger for visibility (scale by cross-section * 50)
        radius = max(self.cross_section_spinbox.value() * 50.0, 200.0)  # At least 200 km for visibility
        sphere.SetRadius(radius)
        sphere.SetThetaResolution(20)
        sphere.SetPhiResolution(20)
        
        object_mapper = vtk.vtkPolyDataMapper()
        object_mapper.SetInputConnection(sphere.GetOutputPort())
        
        if hasattr(self, 'object_actor'):
            self.renderer.RemoveActor(self.object_actor)
            
        self.object_actor = vtk.vtkActor()
        self.object_actor.SetMapper(object_mapper)
        
        # Make it bright red and emissive for high visibility
        self.object_actor.GetProperty().SetColor(1.0, 0.1, 0.1)  # Bright red
        self.object_actor.GetProperty().SetAmbient(0.8)  # High ambient lighting
        self.object_actor.GetProperty().SetDiffuse(0.9)  # High diffuse
        self.object_actor.GetProperty().SetSpecular(0.5)  # Some shininess
        self.object_actor.GetProperty().SetOpacity(1.0)  # Fully opaque
        
        self.renderer.AddActor(self.object_actor)
        
        # Add a trail effect (keep last few positions)
        self.create_satellite_trail()
        
        print(f"Created satellite object with radius {radius:.1f} km")

    def create_satellite_trail(self):
        """Create a trail showing recent satellite positions"""
        if not hasattr(self, 'trail_points'):
            self.trail_points = []
            self.max_trail_length = 20  # Keep last 20 positions
            
        # Create trail polyline
        if len(self.trail_points) > 1:
            trail_vtk_points = vtk.vtkPoints()
            trail_lines = vtk.vtkCellArray()
            
            for point in self.trail_points:
                trail_vtk_points.InsertNextPoint(point)
                
            # Create line segments
            for i in range(len(self.trail_points) - 1):
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, i)
                line.GetPointIds().SetId(1, i + 1)
                trail_lines.InsertNextCell(line)
                
            trail_polydata = vtk.vtkPolyData()
            trail_polydata.SetPoints(trail_vtk_points)
            trail_polydata.SetLines(trail_lines)
            
            trail_mapper = vtk.vtkPolyDataMapper()
            trail_mapper.SetInputData(trail_polydata)
            
            if hasattr(self, 'trail_actor'):
                self.renderer.RemoveActor(self.trail_actor)
                
            self.trail_actor = vtk.vtkActor()
            self.trail_actor.SetMapper(trail_mapper)
            self.trail_actor.GetProperty().SetColor(1.0, 0.5, 0.0)  # Orange trail
            self.trail_actor.GetProperty().SetLineWidth(3.0)
            self.trail_actor.GetProperty().SetOpacity(0.8)
            
            self.renderer.AddActor(self.trail_actor)

    def update_visualization(self):
        """Update the 3D visualization with better satellite tracking"""
        if not self.orbital_path or self.current_time_index >= len(self.orbital_path):
            return
            
        current_point = self.orbital_path[self.current_time_index]
        
        # Update object position
        if hasattr(self, 'object_actor'):
            self.object_actor.SetPosition(current_point.x, current_point.y, current_point.z)
            
            # Update trail
            if not hasattr(self, 'trail_points'):
                self.trail_points = []
                
            self.trail_points.append([current_point.x, current_point.y, current_point.z])
            
            # Limit trail length
            if len(self.trail_points) > self.max_trail_length:
                self.trail_points.pop(0)
                
            # Update trail visualization
            if len(self.trail_points) > 1:
                self.create_satellite_trail()
            
        # Update time slider (but don't trigger callback)
        self.time_slider.blockSignals(True)
        self.time_slider.setValue(self.current_time_index)
        self.time_slider.blockSignals(False)
        
        # Update time label
        hours = current_point.time
        h = int(hours)
        m = int((hours - h) * 60)
        s = int(((hours - h) * 60 - m) * 60)
        self.time_label.setText(f"{h:02d}:{m:02d}:{s:02d}")
        
        # Calculate and display flux
        flux = self.flux_analyzer.analyze_flux_at_point(current_point)
        
        # Update status with more detailed info
        distance_from_earth = np.sqrt(current_point.x**2 + current_point.y**2 + current_point.z**2)
        altitude = distance_from_earth - 6371  # Earth radius
        
        self.status_label.setText(
            f"Time: {current_point.time:.2f}h | Alt: {altitude:.1f}km\n"
            f"Pos: ({current_point.x:.0f}, {current_point.y:.0f}, {current_point.z:.0f}) km\n"
            f"Flux: {flux:.2e} particles/s"
        )
        
        # Force render
        self.vtk_widget.GetRenderWindow().Render()
        
        # Debug print for animation tracking
        if self.current_time_index % 10 == 0:  # Print every 10th step
            print(f"Animation step {self.current_time_index}: t={current_point.time:.2f}h, "
                  f"pos=({current_point.x:.0f},{current_point.y:.0f},{current_point.z:.0f})")

    def animation_step(self):
        """Perform one animation step with wraparound"""
        if not self.orbital_path:
            return
            
        if self.current_time_index >= len(self.orbital_path) - 1:
            # Loop back to beginning for continuous animation
            self.current_time_index = 0
            # Clear trail for new orbit
            if hasattr(self, 'trail_points'):
                self.trail_points.clear()
            print("Animation completed one orbit, restarting...")
        else:
            self.current_time_index += 1
            
        self.update_visualization()
        self.update_plots()

    def start_animation(self):
        """Start the orbital animation with improved feedback"""
        if not self.orbital_path:
            QMessageBox.warning(self, "Warning", "Please load orbital path data first.")
            return
            
        if not self.vtk_data:
            QMessageBox.warning(self, "Warning", "Please load VTK flux data first.")
            return
            
        print(f"Starting animation with {len(self.orbital_path)} orbital points")
        print(f"Current time index: {self.current_time_index}")
        print(f"Animation speed: {self.speed_spinbox.value()} ms")
        
        self.is_playing = True
        self.animation_timer.start(self.speed_spinbox.value())
        
        # Show plot windows if they exist
        if self.slice_window:
            self.slice_window.show()
        if self.spectrum_window:
            self.spectrum_window.show()
        if self.flux_time_window:
            self.flux_time_window.show()
            
        self.play_button.setEnabled(False)
        self.pause_button.setEnabled(True)
        self.stop_button.setEnabled(True)
        
        # Force an immediate update
        self.update_visualization()

    def reset_animation(self):
        """Reset animation to beginning with trail cleanup"""
        self.current_time_index = 0
        
        # Clear trail
        if hasattr(self, 'trail_points'):
            self.trail_points.clear()
            
        # Remove trail actor
        if hasattr(self, 'trail_actor'):
            self.renderer.RemoveActor(self.trail_actor)
            delattr(self, 'trail_actor')
            
        # Clear flux time data
        if self.flux_time_window:
            self.flux_time_window.clear_data()
            
        self.update_visualization()
        print("Animation reset to beginning")

    def jump_to_time(self, time_index):
        """Jump to specific time index with trail update"""
        if 0 <= time_index < len(self.orbital_path):
            old_index = self.current_time_index
            self.current_time_index = time_index
            
            # If jumping backwards or far forward, reset trail
            if time_index < old_index or abs(time_index - old_index) > 5:
                if hasattr(self, 'trail_points'):
                    self.trail_points.clear()
                print(f"Jumped to time index {time_index}, trail reset")
            
            self.update_visualization()
            self.update_plots()

    def update_cross_section(self, radius):
        """Update object cross section and visual representation"""
        self.flux_analyzer.set_cross_section(radius)
        
        # Recreate object with new size
        if hasattr(self, 'object_actor'):
            self.renderer.RemoveActor(self.object_actor)
            self.create_object_representation()
            
            # Update position if we have one
            if self.orbital_path and 0 <= self.current_time_index < len(self.orbital_path):
                current_point = self.orbital_path[self.current_time_index]
                self.object_actor.SetPosition(current_point.x, current_point.y, current_point.z)
                
            self.vtk_widget.GetRenderWindow().Render()
            
        print(f"Updated cross-section radius to {radius:.1f} m")

    def create_path_visualization(self):
        """Create 3D orbital path visualization with improved styling"""
        if not self.orbital_path:
            return
            
        print(f"Creating path visualization with {len(self.orbital_path)} points")
            
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
        self.path_actor.GetProperty().SetColor(0.9, 0.9, 0.2)  # Bright yellow
        self.path_actor.GetProperty().SetLineWidth(3.0)  # Thicker line
        self.path_actor.GetProperty().SetOpacity(0.8)
        
        self.renderer.AddActor(self.path_actor)
        
        # Create object representation
        self.create_object_representation()
        
        # Set initial position
        if self.orbital_path:
            first_point = self.orbital_path[0]
            if hasattr(self, 'object_actor'):
                self.object_actor.SetPosition(first_point.x, first_point.y, first_point.z)
        
        self.vtk_widget.GetRenderWindow().Render()
        print("Path visualization created successfully")
        
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
