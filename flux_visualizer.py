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


    def debug_renderer_state(self):
        """Debug method to check renderer state"""
        print("\n=== RENDERER DEBUG INFO ===")
        
        # Check actors
        actors = self.renderer.GetActors()
        actors.InitTraversal()
        actor_count = 0
        while actors.GetNextActor():
            actor_count += 1
        print(f"Total actors: {actor_count}")
        
        # Check volumes
        volumes = self.renderer.GetVolumes()
        volumes.InitTraversal()
        volume_count = 0
        while volumes.GetNextVolume():
            volume_count += 1
        print(f"Total volumes: {volume_count}")
        
        # Check view props
        props = self.renderer.GetViewProps()
        props.InitTraversal()
        prop_count = 0
        while props.GetNextProp():
            prop_count += 1
        print(f"Total view props: {prop_count}")
        
        # Camera info
        camera = self.renderer.GetActiveCamera()
        print(f"Camera position: {camera.GetPosition()}")
        print(f"Camera focal point: {camera.GetFocalPoint()}")
        print(f"Camera view up: {camera.GetViewUp()}")
        
        # Data bounds if available
        if hasattr(self, 'vtk_data') and self.vtk_data:
            bounds = self.vtk_data.GetBounds()
            print(f"Data bounds: X({bounds[0]:.1f}, {bounds[1]:.1f}) Y({bounds[2]:.1f}, {bounds[3]:.1f}) Z({bounds[4]:.1f}, {bounds[5]:.1f})")
        
        print("=========================\n")

    def test_simple_visualization(self):
        """Test method to create a simple visible object"""
        print("=== TESTING SIMPLE VISUALIZATION ===")
        
        if not self.vtk_data:
            print("No VTK data available")
            return
        
        try:
            # Clear everything first
            #self.clear_field_visualization()
            
            # Create a simple sphere at the center of the data
            bounds = self.vtk_data.GetBounds()
            center = [(bounds[1]+bounds[0])/2, (bounds[3]+bounds[2])/2, (bounds[5]+bounds[4])/2]
            radius = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]) / 10
            
            print(f"Creating test sphere at {center} with radius {radius}")
            
            sphere = vtk.vtkSphereSource()
            sphere.SetCenter(center)
            sphere.SetRadius(radius)
            sphere.SetThetaResolution(20)
            sphere.SetPhiResolution(20)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere.GetOutputPort())
            
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(1.0, 0.0, 0.0)  # Bright red
            actor.GetProperty().SetOpacity(1.0)
            
            self.renderer.AddActor(actor)
            
            # Store as field actor for cleanup
            self.field_actor = actor
            
            # Reset camera and render
            self.renderer.ResetCamera()
            self.vtk_widget.GetRenderWindow().Render()
            
            print("Test sphere created and rendered")
            
        except Exception as e:
            print(f"Error in test visualization: {e}")
            import traceback
            traceback.print_exc()
        
    def setup_ui(self):
        """Setup the user interface with wider control panel"""
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

        self.create_earth_opacity_control(vtk_layout)
        
        # Right panel - controls (WIDER)
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

        debug_group = QGroupBox("Debug Controls")
        debug_layout = QVBoxLayout(debug_group)
        
        self.debug_button = QPushButton("Debug Renderer State")
        self.debug_button.clicked.connect(self.debug_renderer_state)
        debug_layout.addWidget(self.debug_button)
        
        self.test_viz_button = QPushButton("Test Simple Visualization")
        self.test_viz_button.clicked.connect(self.test_simple_visualization)
        debug_layout.addWidget(self.test_viz_button)
        
        # Status
        self.status_label = QLabel("Ready - Load VTK data and orbital CSV to begin")
        self.status_label.setWordWrap(True)
        self.control_layout.addWidget(self.status_label)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.control_layout.addWidget(self.progress_bar)
        
        self.control_layout.addStretch()
        
        # WIDER control panel (was 350, now 450)
        control_panel.setMaximumWidth(450)
        control_panel.setMinimumWidth(450)
        splitter.addWidget(control_panel)
        
        # Set splitter proportions (adjusted for wider panel)
        splitter.setSizes([750, 450])
        
        self.connect_signals()
        
    def setup_vtk(self):
        """Setup VTK 3D visualization - DEBUG VERSION"""
        print("=== SETTING UP VTK RENDERER ===")
        
        # Renderer
        self.renderer = vtk.vtkRenderer()
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.renderer.SetBackground(0.05, 0.05, 0.15)  # Dark blue background
        
        # Camera setup for Earth-centered view
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(20000, 20000, 10000)  # km
        camera.SetFocalPoint(0, 0, 0)
        camera.SetViewUp(0, 0, 1)
        
        print("Renderer and camera setup complete")
        
        # Create Earth representation ONCE, with error handling
        print(">>> About to call create_earth_representation...")
        earth_success = self.create_earth_representation()
        print(f">>> create_earth_representation returned: {earth_success} (type: {type(earth_success)})")
        
        if earth_success:
            print("Earth created successfully")
        else:
            print("WARNING: Earth creation failed, continuing without Earth")
        
        # Add coordinate axes
        print("Setting up coordinate axes...")
        try:
            self.setup_coordinate_axes()
            print("Coordinate axes setup complete")
        except Exception as e:
            print(f"WARNING: Coordinate axes setup failed: {e}")
        
        print("=== VTK SETUP COMPLETE ===")

    def setup_coordinate_axes(self):
        """Setup dynamic coordinate axes with snap-to-axis controls"""
        # Create axes actor with proper labels and units
        self.axes_actor = vtk.vtkAxesActor()
        
        # Set axis labels with units
        self.axes_actor.SetXAxisLabelText("X (km)")
        self.axes_actor.SetYAxisLabelText("Y (km)")
        self.axes_actor.SetZAxisLabelText("Z (km)")
        
        # Style the axes
        self.axes_actor.SetTotalLength(5000, 5000, 5000)  # 5000 km length
        self.axes_actor.SetShaftTypeToLine()
        self.axes_actor.SetAxisLabels(True)
        
        # Improve text properties
        for axis_label in [self.axes_actor.GetXAxisCaptionActor2D(),
                          self.axes_actor.GetYAxisCaptionActor2D(),
                          self.axes_actor.GetZAxisCaptionActor2D()]:
            axis_label.GetTextActor().SetTextScaleModeToNone()
            axis_label.GetCaptionTextProperty().SetFontSize(12)
            axis_label.GetCaptionTextProperty().SetColor(1, 1, 1)
            axis_label.GetCaptionTextProperty().BoldOn()
        
        # Create orientation widget to show axes in corner
        self.orientation_widget = vtk.vtkOrientationMarkerWidget()
        self.orientation_widget.SetOrientationMarker(self.axes_actor)
        self.orientation_widget.SetInteractor(self.vtk_widget)
        self.orientation_widget.SetViewport(0.0, 0.0, 0.3, 0.3)  # Bottom left corner
        self.orientation_widget.EnabledOn()
        self.orientation_widget.InteractiveOff()  # Don't allow user to move it
        
        # Add snap-to-axis buttons next to the orientation widget
        self.setup_snap_to_axis_buttons()
        
        print("Dynamic coordinate axes with snap buttons added")

    def setup_snap_to_axis_buttons(self):
        """Add snap-to-axis buttons overlaid on the VTK widget"""
        # Create a widget to hold the snap buttons
        self.snap_buttons_widget = QWidget(self.vtk_widget)
        self.snap_buttons_widget.setGeometry(10, 10, 200, 30)  # Top-left corner
        
        # Create horizontal layout for buttons
        snap_layout = QHBoxLayout(self.snap_buttons_widget)
        snap_layout.setContentsMargins(2, 2, 2, 2)
        snap_layout.setSpacing(3)
        
        # Create snap buttons with tooltips
        self.snap_x_button = QPushButton("X")
        self.snap_x_button.setMaximumSize(25, 25)
        self.snap_x_button.setToolTip("Snap to X-axis view (X→right, Y→up, Z→out)")
        self.snap_x_button.clicked.connect(self.snap_to_x_axis)
        self.snap_x_button.setStyleSheet("""
            QPushButton {
                background-color: rgba(200, 50, 50, 180);
                color: white;
                border: 1px solid white;
                border-radius: 3px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: rgba(255, 70, 70, 200);
            }
        """)
        
        self.snap_y_button = QPushButton("Y")
        self.snap_y_button.setMaximumSize(25, 25)
        self.snap_y_button.setToolTip("Snap to Y-axis view (Y→right, Z→up, X→out)")
        self.snap_y_button.clicked.connect(self.snap_to_y_axis)
        self.snap_y_button.setStyleSheet("""
            QPushButton {
                background-color: rgba(50, 200, 50, 180);
                color: white;
                border: 1px solid white;
                border-radius: 3px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: rgba(70, 255, 70, 200);
            }
        """)
        
        self.snap_z_button = QPushButton("Z")
        self.snap_z_button.setMaximumSize(25, 25)
        self.snap_z_button.setToolTip("Snap to Z-axis view (X→right, Z→up, Y→out)")
        self.snap_z_button.clicked.connect(self.snap_to_z_axis)
        self.snap_z_button.setStyleSheet("""
            QPushButton {
                background-color: rgba(50, 50, 200, 180);
                color: white;
                border: 1px solid white;
                border-radius: 3px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: rgba(70, 70, 255, 200);
            }
        """)
        
        # Add ISO button for isometric view
        self.snap_iso_button = QPushButton("ISO")
        self.snap_iso_button.setMaximumSize(35, 25)
        self.snap_iso_button.setToolTip("Snap to isometric 3D view")
        self.snap_iso_button.clicked.connect(self.snap_to_isometric)
        self.snap_iso_button.setStyleSheet("""
            QPushButton {
                background-color: rgba(100, 100, 100, 180);
                color: white;
                border: 1px solid white;
                border-radius: 3px;
                font-weight: bold;
                font-size: 9px;
            }
            QPushButton:hover {
                background-color: rgba(150, 150, 150, 200);
            }
        """)
        
        # Add buttons to layout
        snap_layout.addWidget(QLabel("View:"))
        snap_layout.addWidget(self.snap_x_button)
        snap_layout.addWidget(self.snap_y_button) 
        snap_layout.addWidget(self.snap_z_button)
        snap_layout.addWidget(self.snap_iso_button)
        snap_layout.addStretch()
        
        # Make the widget visible
        self.snap_buttons_widget.show()
        self.snap_buttons_widget.raise_()

    def snap_to_x_axis(self):
        """Snap camera to X-axis aligned view"""
        print("Snapping to X-axis view...")
        
        if not hasattr(self, 'vtk_data') or not self.vtk_data:
            # Use default center
            focal_point = [0, 0, 0]
            distance = 30000
        else:
            # Use data center
            bounds = self.vtk_data.GetBounds()
            focal_point = [
                (bounds[1] + bounds[0]) / 2,
                (bounds[3] + bounds[2]) / 2,
                (bounds[5] + bounds[4]) / 2
            ]
            max_range = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
            distance = max_range * 1.5
        
        camera = self.renderer.GetActiveCamera()
        camera.SetFocalPoint(focal_point)
        
        # X→right, Y→up, Z→out of page
        # Camera looks from negative Z direction
        camera.SetPosition(focal_point[0], focal_point[1], focal_point[2] - distance)
        camera.SetViewUp(0, 1, 0)  # Y axis points up
        
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()

    def snap_to_y_axis(self):
        """Snap camera to Y-axis aligned view"""
        print("Snapping to Y-axis view...")
        
        if not hasattr(self, 'vtk_data') or not self.vtk_data:
            focal_point = [0, 0, 0]
            distance = 30000
        else:
            bounds = self.vtk_data.GetBounds()
            focal_point = [
                (bounds[1] + bounds[0]) / 2,
                (bounds[3] + bounds[2]) / 2,
                (bounds[5] + bounds[4]) / 2
            ]
            max_range = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
            distance = max_range * 1.5
        
        camera = self.renderer.GetActiveCamera()
        camera.SetFocalPoint(focal_point)
        
        # Y→right, Z→up, X→out of page
        # Camera looks from negative X direction
        camera.SetPosition(focal_point[0] - distance, focal_point[1], focal_point[2])
        camera.SetViewUp(0, 0, 1)  # Z axis points up
        
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()

    def snap_to_z_axis(self):
        """Snap camera to Z-axis aligned view"""
        print("Snapping to Z-axis view...")
        
        if not hasattr(self, 'vtk_data') or not self.vtk_data:
            focal_point = [0, 0, 0]
            distance = 30000
        else:
            bounds = self.vtk_data.GetBounds()
            focal_point = [
                (bounds[1] + bounds[0]) / 2,
                (bounds[3] + bounds[2]) / 2,
                (bounds[5] + bounds[4]) / 2
            ]
            max_range = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
            distance = max_range * 1.5
        
        camera = self.renderer.GetActiveCamera()
        camera.SetFocalPoint(focal_point)
        
        # X→right, Z→up, Y→out of page
        # Camera looks from negative Y direction
        camera.SetPosition(focal_point[0], focal_point[1] - distance, focal_point[2])
        camera.SetViewUp(0, 0, 1)  # Z axis points up
        
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()

    def snap_to_isometric(self):
        """Snap camera to nice isometric 3D view"""
        print("Snapping to isometric view...")
        
        if not hasattr(self, 'vtk_data') or not self.vtk_data:
            focal_point = [0, 0, 0]
            distance = 30000
        else:
            bounds = self.vtk_data.GetBounds()
            focal_point = [
                (bounds[1] + bounds[0]) / 2,
                (bounds[3] + bounds[2]) / 2,
                (bounds[5] + bounds[4]) / 2
            ]
            max_range = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
            distance = max_range * 1.5
        
        camera = self.renderer.GetActiveCamera()
        camera.SetFocalPoint(focal_point)
        
        # Nice isometric view: roughly 45° angles
        camera.SetPosition(
            focal_point[0] + distance * 0.7,
            focal_point[1] + distance * 0.7,
            focal_point[2] + distance * 0.5
        )
        camera.SetViewUp(0, 0, 1)  # Z axis points up
        
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()

    def create_earth_representation(self):
        """Create Earth with texture, lat/long grid, and opacity control"""
        print("Creating Earth with texture, grid, and controls...")
        
        try:
            # Remove existing Earth actors if they exist
            self.cleanup_existing_earth_actors()
            
            # 1. Create Earth sphere
            print("Creating Earth sphere...")
            earth_sphere = vtk.vtkSphereSource()
            earth_sphere.SetRadius(6371.0)  # Earth radius in km
            earth_sphere.SetThetaResolution(120)  # Higher resolution to reduce artifacts
            earth_sphere.SetPhiResolution(120)
            earth_sphere.SetCenter(0.0, 0.0, 0.0)
            earth_sphere.Update()
            
            # 2. Manually add correct texture coordinates
            sphere_data = earth_sphere.GetOutput()
            self.add_correct_texture_coordinates(sphere_data)
            
            # 3. Create mapper directly from sphere data
            earth_mapper = vtk.vtkPolyDataMapper()
            earth_mapper.SetInputData(sphere_data)
            
            # 4. Create Earth actor
            self.earth_actor = vtk.vtkActor()
            self.earth_actor.SetMapper(earth_mapper)
            
            # 5. Try to load and apply Earth texture
            earth_texture = self.load_earth_texture()
            if earth_texture:
                earth_texture.SetRepeat(0)  # Don't repeat texture
                earth_texture.SetInterpolate(1)  # Smooth interpolation
                # Additional settings to reduce artifacts
                earth_texture.SetWrap(vtk.vtkTexture.ClampToEdge)
                self.earth_actor.SetTexture(earth_texture)
                print("Earth texture applied with artifact fixes")
            else:
                # Fallback to solid color
                print("Using fallback solid Earth color")
                self.earth_actor.GetProperty().SetColor(0.2, 0.5, 0.8)
            
            # 6. Set Earth properties with default opacity
            earth_property = self.earth_actor.GetProperty()
            earth_property.SetOpacity(0.8)  # Default 80% opacity
            earth_property.SetAmbient(0.4)
            earth_property.SetDiffuse(0.8)
            earth_property.SetSpecular(0.05)
            earth_property.SetSpecularPower(10)
            
            self.earth_actor.SetVisibility(True)
            self.renderer.AddActor(self.earth_actor)
            
            # 7. Create latitude lines (every 15 degrees for better coverage)
            print("Creating latitude grid lines...")
            self.lat_long_actors = []
            
            for lat_deg in range(-75, 90, 15):  # Every 15 degrees, -75 to 75
                lat_rad = np.radians(lat_deg)
                radius = 6372.0 * np.cos(lat_rad)  # Slightly above Earth surface
                height = 6372.0 * np.sin(lat_rad)
                
                if radius > 100:  # Skip very small circles near poles
                    circle = self.create_circle(radius, height, 'lat')
                    if circle:
                        self.lat_long_actors.append(circle)
                        self.renderer.AddActor(circle)
            
            # 8. Create longitude lines (every 15 degrees for full coverage)
            print("Creating longitude grid lines...")
            for lon_deg in range(0, 360, 15):  # Every 15 degrees, full 360°
                meridian = self.create_meridian(lon_deg)
                if meridian:
                    self.lat_long_actors.append(meridian)
                    self.renderer.AddActor(meridian)
            
            # 9. Create highlighted equator and prime meridian
            print("Creating equator and prime meridian...")
            self.equator_actor = self.create_circle(6372.5, 0.0, 'equator')
            if self.equator_actor:
                self.renderer.AddActor(self.equator_actor)
            
            # Create prime meridian (highlighted)
            self.prime_meridian_actor = self.create_meridian(0, highlight=True)
            if self.prime_meridian_actor:
                self.renderer.AddActor(self.prime_meridian_actor)
            
            # 10. Add Earth opacity control to UI
            #self.add_earth_opacity_control()
            
            # Force render
            if hasattr(self, 'vtk_widget'):
                self.vtk_widget.GetRenderWindow().Render()
            
            print("Earth with texture, grid, and opacity control created successfully")
            return True
            
        except Exception as e:
            print(f"Error creating enhanced Earth: {e}")
            import traceback
            traceback.print_exc()
            return False

    def create_earth_opacity_control(self, parent_layout):
        """Create Earth opacity control and add it to the parent layout"""
        try:
            # Create opacity control widget
            opacity_widget = QWidget()
            opacity_layout = QHBoxLayout(opacity_widget)
            opacity_layout.setContentsMargins(10, 5, 10, 5)
            opacity_layout.setSpacing(10)
            
            # Add label
            opacity_label = QLabel("Earth Opacity:")
            opacity_label.setStyleSheet("color: white; font-weight: bold;")
            opacity_layout.addWidget(opacity_label)
            
            # Add slider
            self.earth_opacity_slider = QSlider(Qt.Orientation.Horizontal)
            self.earth_opacity_slider.setRange(0, 100)
            self.earth_opacity_slider.setValue(80)  # Default 80%
            self.earth_opacity_slider.setFixedWidth(150)
            self.earth_opacity_slider.valueChanged.connect(self.update_earth_opacity)
            opacity_layout.addWidget(self.earth_opacity_slider)
            
            # Add value label
            self.earth_opacity_value_label = QLabel("80%")
            self.earth_opacity_value_label.setStyleSheet("color: white; font-weight: bold; min-width: 35px;")
            opacity_layout.addWidget(self.earth_opacity_value_label)
            
            # Add spacer
            opacity_layout.addStretch()
            
            # Style the widget
            opacity_widget.setStyleSheet("""
                QWidget {
                    background-color: rgba(40, 40, 40, 180);
                    border: 1px solid #555;
                    border-radius: 5px;
                    margin: 2px;
                }
                QSlider::groove:horizontal {
                    border: 1px solid #666;
                    height: 6px;
                    background: #333;
                    border-radius: 3px;
                }
                QSlider::handle:horizontal {
                    background: #888;
                    border: 1px solid #555;
                    width: 16px;
                    margin: -5px 0;
                    border-radius: 8px;
                }
                QSlider::handle:horizontal:hover {
                    background: #aaa;
                }
            """)
            
            opacity_widget.setFixedHeight(35)
            
            # Add to parent layout
            parent_layout.addWidget(opacity_widget)
            
            print("Earth opacity control created successfully")
            
        except Exception as e:
            print(f"Error creating Earth opacity control: {e}")

    def update_earth_opacity(self, value):
        """Update Earth opacity from slider"""
        try:
            opacity = value / 100.0
            self.earth_opacity_value_label.setText(f"{value}%")
            
            if hasattr(self, 'earth_actor') and self.earth_actor:
                self.earth_actor.GetProperty().SetOpacity(opacity)
                
                # Also update grid line opacity proportionally
                grid_opacity = min(0.8, opacity + 0.2)  # Keep grid slightly more visible
                
                if hasattr(self, 'lat_long_actors'):
                    for actor in self.lat_long_actors:
                        if actor:
                            actor.GetProperty().SetOpacity(grid_opacity * 0.6)
                
                if hasattr(self, 'equator_actor') and self.equator_actor:
                    self.equator_actor.GetProperty().SetOpacity(grid_opacity * 0.8)
                    
                if hasattr(self, 'prime_meridian_actor') and self.prime_meridian_actor:
                    self.prime_meridian_actor.GetProperty().SetOpacity(grid_opacity * 0.8)
                
                # Force render
                if hasattr(self, 'vtk_widget'):
                    self.vtk_widget.GetRenderWindow().Render()
                    
        except Exception as e:
            print(f"Error updating Earth opacity: {e}")
        
    def create_meridian(self, longitude_deg, highlight=False):
        """Create a meridian (longitude line) from pole to pole"""
        try:
            lon_rad = np.radians(longitude_deg)
            
            # Create points from south pole to north pole
            n_points = 120  # Higher resolution for smoother lines
            latitudes = np.linspace(-np.pi/2, np.pi/2, n_points)
            
            points = vtk.vtkPoints()
            lines = vtk.vtkCellArray()
            
            for i, lat in enumerate(latitudes):
                x = 6372.0 * np.cos(lat) * np.cos(lon_rad)  # 1km above Earth
                y = 6372.0 * np.cos(lat) * np.sin(lon_rad)
                z = 6372.0 * np.sin(lat)
                points.InsertNextPoint(x, y, z)
                
                if i > 0:
                    line = vtk.vtkLine()
                    line.GetPointIds().SetId(0, i-1)
                    line.GetPointIds().SetId(1, i)
                    lines.InsertNextCell(line)
            
            # Create polydata
            meridian_polydata = vtk.vtkPolyData()
            meridian_polydata.SetPoints(points)
            meridian_polydata.SetLines(lines)
            
            # Create mapper
            meridian_mapper = vtk.vtkPolyDataMapper()
            meridian_mapper.SetInputData(meridian_polydata)
            
            # Create actor
            meridian_actor = vtk.vtkActor()
            meridian_actor.SetMapper(meridian_mapper)
            
            # Set properties based on whether it's highlighted
            if highlight:  # Prime meridian or special meridian
                meridian_actor.GetProperty().SetColor(1.0, 0.6, 0.0)  # Orange
                meridian_actor.GetProperty().SetLineWidth(2.5)
                meridian_actor.GetProperty().SetOpacity(0.8)
            else:  # Regular meridian
                meridian_actor.GetProperty().SetColor(0.7, 0.7, 0.7)  # Light gray
                meridian_actor.GetProperty().SetLineWidth(1.5)
                meridian_actor.GetProperty().SetOpacity(0.6)
            
            return meridian_actor
            
        except Exception as e:
            print(f"Error creating meridian at {longitude_deg}°: {e}")
            return None

    def add_correct_texture_coordinates(self, sphere_data):
        """Add texture coordinates with proper seam handling to eliminate artifacts"""
        try:
            print("Computing texture coordinates with seam artifact fix...")
            
            points = sphere_data.GetPoints()
            num_points = points.GetNumberOfPoints()
            
            # Create texture coordinate array
            tex_coords = vtk.vtkFloatArray()
            tex_coords.SetNumberOfComponents(2)
            tex_coords.SetNumberOfTuples(num_points)
            tex_coords.SetName("TextureCoordinates")
            
            for i in range(num_points):
                # Get 3D point
                point = points.GetPoint(i)
                x, y, z = point
                
                # Convert to spherical coordinates
                r = np.sqrt(x*x + y*y + z*z)
                
                # Calculate longitude (u coordinate) with special seam handling
                longitude = np.arctan2(y, x)  # -π to π
                
                # Handle the seam at ±π (International Date Line)
                # This is the key fix for the Pacific artifact
                if longitude < 0:
                    longitude += 2 * np.pi  # Convert to 0 to 2π
                    
                u = longitude / (2 * np.pi)  # 0 to 1
                
                # Add small epsilon to prevent exact 0 or 1 values that cause seams
                epsilon = 1e-6
                u = max(epsilon, min(1.0 - epsilon, u))
                
                # Calculate latitude (v coordinate)
                latitude = np.arcsin(np.clip(z / r, -1.0, 1.0))  # -π/2 to π/2
                
                # Map latitude to texture coordinates (0 to 1)
                v = (latitude + np.pi/2) / np.pi  # 0 to 1
                v = max(epsilon, min(1.0 - epsilon, v))
                
                # Set texture coordinates
                tex_coords.SetTuple2(i, u, v)
            
            # Add texture coordinates to the sphere data
            sphere_data.GetPointData().SetTCoords(tex_coords)
            print(f"Added seam-fixed texture coordinates for {num_points} points")
            
        except Exception as e:
            print(f"Error computing texture coordinates: {e}")
            import traceback
            traceback.print_exc()

    def load_earth_texture(self):
        """Load Earth texture with proper handling"""
        try:
            # List of possible texture files
            texture_files = [
                "earth_texture.jpg", "earth_map.jpg", "world_map.jpg",
                "earth_texture.png", "earth_map.png", "world_map.png",
                "blue_marble.jpg", "blue_marble.png",
                "natural_earth.jpg", "natural_earth.png",
                "earth.jpg", "earth.png", "world.jpg", "world.png"
            ]
            
            # Try to find and load texture file
            for filename in texture_files:
                texture = self.try_load_texture_file(filename)
                if texture:
                    print(f"Successfully loaded texture: {filename}")
                    return texture
            
            # If no local file found, create procedural texture
            print("No texture file found, creating procedural Earth texture...")
            return self.create_procedural_earth_texture()
            
        except Exception as e:
            print(f"Error loading Earth texture: {e}")
            return None

    def try_load_texture_file(self, filename):
        """Try to load a specific texture file with settings to prevent artifacts"""
        try:
            import os
            
            if not os.path.exists(filename):
                return None
            
            # Create appropriate reader
            if filename.lower().endswith(('.jpg', '.jpeg')):
                reader = vtk.vtkJPEGReader()
            elif filename.lower().endswith('.png'):
                reader = vtk.vtkPNGReader()
            else:
                return None
            
            reader.SetFileName(filename)
            reader.Update()
            
            # Check if image loaded
            if reader.GetOutput().GetNumberOfPoints() == 0:
                return None
            
            # Create texture with specific settings to prevent artifacts
            texture = vtk.vtkTexture()
            texture.SetInputConnection(reader.GetOutputPort())
            texture.InterpolateOn()  # Smooth interpolation
            texture.RepeatOff()      # Don't repeat - this is crucial
            texture.EdgeClampOn()    # Clamp edges
            
            # Additional settings to prevent seam artifacts
            texture.SetWrap(vtk.vtkTexture.ClampToEdge)
            texture.SetQualityTo32Bit()  # Higher quality
            
            print(f"Loaded texture {filename} with anti-artifact settings")
            return texture
            
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            return None

    def create_procedural_earth_texture(self):
        """Create a simple procedural Earth texture (improved)"""
        try:
            print("Creating improved procedural Earth texture...")
            
            # Standard equirectangular dimensions
            width, height = 1024, 512
            
            # Create image data
            image = vtk.vtkImageData()
            image.SetDimensions(width, height, 1)
            image.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 3)  # RGB
            
            # Create a simple but recognizable Earth pattern
            for y in range(height):
                for x in range(width):
                    # Convert pixel to lat/lon
                    lon = (x / width) * 360 - 180  # -180 to 180
                    lat = 90 - (y / height) * 180   # 90 to -90
                    
                    # Simple continent approximation
                    is_land = self.simple_land_check(lon, lat)
                    
                    if is_land:
                        # Land colors (green/brown with variation)
                        variation = (x + y) % 50 / 50.0
                        r = int(80 + variation * 60)   # Brown-green
                        g = int(100 + variation * 80)  # Green
                        b = int(40 + variation * 40)   # Minimal blue
                    else:
                        # Ocean colors (blue with variation)
                        variation = (x * 3 + y * 2) % 40 / 40.0
                        r = int(20 + variation * 30)   # Minimal red
                        g = int(60 + variation * 60)   # Medium green
                        b = int(120 + variation * 80)  # Strong blue
                    
                    # Clamp values
                    r = max(0, min(255, r))
                    g = max(0, min(255, g))
                    b = max(0, min(255, b))
                    
                    # Set pixel
                    image.SetScalarComponentFromFloat(x, y, 0, 0, r)
                    image.SetScalarComponentFromFloat(x, y, 0, 1, g)
                    image.SetScalarComponentFromFloat(x, y, 0, 2, b)
            
            # Create texture
            texture = vtk.vtkTexture()
            texture.SetInputData(image)
            texture.InterpolateOn()
            texture.RepeatOff()
            texture.EdgeClampOn()
            
            print("Procedural Earth texture created")
            return texture
            
        except Exception as e:
            print(f"Error creating procedural texture: {e}")
            return None

    def simple_land_check(self, lon, lat):
        """Simple land/ocean check for procedural texture"""
        # Very simplified continent shapes
        # North America
        if (-140 < lon < -60 and 15 < lat < 70):
            return True
        # South America  
        if (-80 < lon < -35 and -55 < lat < 15):
            return True
        # Africa
        if (-20 < lon < 50 and -35 < lat < 35):
            return True
        # Europe
        if (-10 < lon < 40 and 35 < lat < 70):
            return True
        # Asia
        if (40 < lon < 180 and 10 < lat < 75):
            return True
        # Australia
        if (110 < lon < 155 and -45 < lat < -10):
            return True
        
        return False

    # def add_earth_opacity_control(self):
    #     """Add Earth opacity control below the visualization"""
    #     try:
    #         # Find the main layout that contains the VTK widget
    #         vtk_frame = self.vtk_widget.parent()
    #         vtk_layout = vtk_frame.layout()
            
    #         # Create opacity control widget
    #         opacity_control = QWidget()
    #         opacity_layout = QHBoxLayout(opacity_control)
    #         opacity_layout.setContentsMargins(10, 5, 10, 5)
            
    #         # Add label
    #         opacity_label = QLabel("Earth Opacity:")
    #         opacity_layout.addWidget(opacity_label)
            
    #         # Add slider
    #         self.earth_opacity_slider = QSlider(Qt.Orientation.Horizontal)
    #         self.earth_opacity_slider.setRange(0, 100)
    #         self.earth_opacity_slider.setValue(80)  # Default 80%
    #         self.earth_opacity_slider.setMaximumWidth(200)
    #         self.earth_opacity_slider.valueChanged.connect(self.update_earth_opacity)
    #         opacity_layout.addWidget(self.earth_opacity_slider)
            
    #         # Add value label
    #         self.earth_opacity_label = QLabel("80%")
    #         self.earth_opacity_label.setMinimumWidth(40)
    #         opacity_layout.addWidget(self.earth_opacity_label)
            
    #         # Add some space
    #         opacity_layout.addStretch()
            
    #         # Style the control
    #         opacity_control.setStyleSheet("""
    #             QWidget {
    #                 background-color: rgba(50, 50, 50, 200);
    #                 border: 1px solid #666;
    #                 border-radius: 5px;
    #             }
    #             QLabel {
    #                 color: white;
    #                 font-weight: bold;
    #             }
    #             QSlider::groove:horizontal {
    #                 border: 1px solid #999;
    #                 height: 8px;
    #                 background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #B0B0B0, stop:1 #A0A0A0);
    #                 margin: 2px 0;
    #                 border-radius: 4px;
    #             }
    #             QSlider::handle:horizontal {
    #                 background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #b4b4b4, stop:1 #8f8f8f);
    #                 border: 1px solid #5c5c5c;
    #                 width: 18px;
    #                 margin: -2px 0;
    #                 border-radius: 3px;
    #             }
    #         """)
            
    #         # Add to the VTK layout at the bottom
    #         vtk_layout.addWidget(opacity_control)
            
    #         print("Earth opacity control added below visualization")
            
    #     except Exception as e:
    #         print(f"Error adding opacity control: {e}")
    #         # Fallback: add to main control panel if VTK layout fails
    #         try:
    #             self.add_opacity_control_to_sidebar()
    #         except:
    #             print("Could not add opacity control anywhere")

    # def add_opacity_control_to_sidebar(self):
    #     """Fallback: add opacity control to sidebar"""
    #     # Find analysis parameters group
    #     for i in range(self.control_layout.count()):
    #         item = self.control_layout.itemAt(i)
    #         if item and isinstance(item.widget(), QGroupBox) and item.widget().title() == "Analysis Parameters":
    #             analysis_group = item.widget()
    #             analysis_layout = analysis_group.layout()
                
    #             # Add Earth opacity control
    #             self.earth_opacity_slider = QSlider(Qt.Orientation.Horizontal)
    #             self.earth_opacity_slider.setRange(0, 100)
    #             self.earth_opacity_slider.setValue(80)
    #             self.earth_opacity_slider.valueChanged.connect(self.update_earth_opacity)
                
    #             self.earth_opacity_label = QLabel("80%")
    #             self.earth_opacity_label.setMinimumWidth(50)
                
    #             analysis_layout.addRow("Earth Opacity:", self.earth_opacity_slider)
    #             analysis_layout.addRow("", self.earth_opacity_label)
    #            break

    # def update_earth_opacity(self, value):
    #     """Update Earth opacity from slider"""
    #     try:
    #         opacity = value / 100.0
    #         self.earth_opacity_label.setText(f"{value}%")
            
    #         if hasattr(self, 'earth_actor') and self.earth_actor:
    #             self.earth_actor.GetProperty().SetOpacity(opacity)
                
    #             # Force render
    #             if hasattr(self, 'vtk_widget'):
    #                 self.vtk_widget.GetRenderWindow().Render()
                    
    #             print(f"Earth opacity updated to {value}%")
                
    #     except Exception as e:
    #         print(f"Error updating Earth opacity: {e}")
    
    def cleanup_existing_earth_actors(self):
        """Clean up all existing Earth-related actors"""
        actors_to_cleanup = [
            'earth_actor', 'ocean_actor', 'earth_wireframe_actor', 
            'equator_actor', 'prime_meridian_actor'
        ]
        
        for actor_name in actors_to_cleanup:
            if hasattr(self, actor_name):
                actor = getattr(self, actor_name)
                if actor:
                    self.renderer.RemoveActor(actor)
                    setattr(self, actor_name, None)
        
        # Clean up actor lists
        if hasattr(self, 'lat_long_actors') and self.lat_long_actors:
            for actor in self.lat_long_actors:
                if actor:
                    self.renderer.RemoveActor(actor)
            self.lat_long_actors = []

    def download_earth_texture(self):
        """Download a free Earth texture (helper function)"""
        print("""
To get a high-quality Earth texture, you can download one from:

1. NASA Blue Marble: 
   https://visibleearth.nasa.gov/images/57752/blue-marble-land-surface-shallow-water-and-shaded-topography
   
2. Natural Earth: 
   https://www.naturalearthdata.com/downloads/10m-raster-data/
   
3. Free texture sites:
   - https://www.solarsystemscope.com/textures/ (Planet textures)
   - https://planetpixelemporium.com/planets.html (Free planet textures)

Save the image as 'earth_texture.jpg' in the same directory as your script.
For best results, use an equirectangular projection (2:1 aspect ratio).
        """)

    def cleanup_existing_earth_actors(self):
        """Clean up all existing Earth-related actors"""
        actors_to_cleanup = [
            'earth_actor', 'ocean_actor', 'earth_wireframe_actor', 
            'equator_actor'
        ]
        
        for actor_name in actors_to_cleanup:
            if hasattr(self, actor_name):
                actor = getattr(self, actor_name)
                if actor:
                    self.renderer.RemoveActor(actor)
                    setattr(self, actor_name, None)
        
        # Clean up actor lists
        if hasattr(self, 'lat_long_actors') and self.lat_long_actors:
            for actor in self.lat_long_actors:
                if actor:
                    self.renderer.RemoveActor(actor)
            self.lat_long_actors = []
        
    def clear_field_visualization(self):
        """Clear field visualization while preserving Earth texture and lat/long grid"""
        print("=== CLEARING FIELD VISUALIZATION (PRESERVING EARTH) ===")
        
        actors_removed = 0
        
        # Store Earth-related actors that should be preserved
        earth_actors_to_preserve = {
            'earth_actor': getattr(self, 'earth_actor', None),
            'equator_actor': getattr(self, 'equator_actor', None),
            'lat_long_actors': getattr(self, 'lat_long_actors', [])
        }
        
        # Flatten the list of actors to preserve for easy checking
        preserve_list = []
        if earth_actors_to_preserve['earth_actor']:
            preserve_list.append(earth_actors_to_preserve['earth_actor'])
        if earth_actors_to_preserve['equator_actor']:
            preserve_list.append(earth_actors_to_preserve['equator_actor'])
        preserve_list.extend(earth_actors_to_preserve['lat_long_actors'])
        
        # Remove field visualization actors (but preserve Earth)
        if hasattr(self, 'field_actor') and self.field_actor:
            if self.field_actor not in preserve_list:
                self.renderer.RemoveActor(self.field_actor)
                self.field_actor = None
                actors_removed += 1
                print("Removed field_actor")
        
        # Remove volume actor
        if hasattr(self, 'volume_actor') and self.volume_actor:
            self.renderer.RemoveVolume(self.volume_actor)
            self.volume_actor = None
            actors_removed += 1
            print("Removed volume_actor")
        
        # Remove other visualization actors while preserving Earth
        visualization_actor_lists = ['wireframe_actors', 'slice_actors', 'isosurface_actors']
        
        for actor_list_name in visualization_actor_lists:
            if hasattr(self, actor_list_name):
                actor_list = getattr(self, actor_list_name)
                for actor in actor_list:
                    if actor and actor not in preserve_list:
                        self.renderer.RemoveActor(actor)
                        actors_removed += 1
                setattr(self, actor_list_name, [])
                if actor_list:
                    print(f"Removed {actor_list_name}")
        
        # Remove scalar bar
        if hasattr(self, 'scalar_bar') and self.scalar_bar:
            self.renderer.RemoveViewProp(self.scalar_bar)
            self.scalar_bar = None
            print("Removed scalar_bar")
        
        # Verify Earth actors are still in renderer (re-add if missing)
        self.verify_earth_actors_in_renderer(earth_actors_to_preserve)
        
        print(f"Total field actors removed: {actors_removed} (Earth preserved)")
        print("=== CLEAR COMPLETE ===")
    
    def verify_earth_actors_in_renderer(self, earth_actors_to_preserve):
        """Ensure all Earth actors are still in the renderer"""
        # Get current actors in renderer
        current_actors = []
        actors = self.renderer.GetActors()
        actors.InitTraversal()
        while True:
            actor = actors.GetNextActor()
            if not actor:
                break
            current_actors.append(actor)
        
        # Check and re-add Earth actor if missing
        if earth_actors_to_preserve['earth_actor']:
            if earth_actors_to_preserve['earth_actor'] not in current_actors:
                self.renderer.AddActor(earth_actors_to_preserve['earth_actor'])
                print("Re-added Earth texture actor")
        
        # Check and re-add equator actor if missing
        if earth_actors_to_preserve['equator_actor']:
            if earth_actors_to_preserve['equator_actor'] not in current_actors:
                self.renderer.AddActor(earth_actors_to_preserve['equator_actor'])
                print("Re-added equator actor")
        
        # Check and re-add lat/long grid actors if missing
        readded_grid_count = 0
        for grid_actor in earth_actors_to_preserve['lat_long_actors']:
            if grid_actor and grid_actor not in current_actors:
                self.renderer.AddActor(grid_actor)
                readded_grid_count += 1
        
        if readded_grid_count > 0:
            print(f"Re-added {readded_grid_count} lat/long grid actors")
        
    def setup_visualization_controls(self):
        """Enhanced controls with descriptions UNDER slider rows"""
        
        # Find the control layout
        control_layout = None
        for i in range(self.control_layout.count()):
            item = self.control_layout.itemAt(i)
            if item and isinstance(item.widget(), QGroupBox) and item.widget().title() == "Analysis Parameters":
                control_layout = self.control_layout
                insert_index = i + 1
                break
        
        if not control_layout:
            control_layout = self.control_layout
            insert_index = -1
        
        # Create visualization controls group
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
        self.viz_mode_combo.setCurrentText("Point Cloud")
        self.viz_mode_combo.currentTextChanged.connect(self.change_visualization_mode)
        mode_layout.addRow("Visualization Mode:", self.viz_mode_combo)
        
        # Wireframe Style Controls (only visible when Wireframe is selected)
        self.wireframe_controls = QWidget()
        wireframe_layout = QVBoxLayout(self.wireframe_controls)  # Changed to VBox for descriptions
        
        # Wireframe style selector
        style_row = QHBoxLayout()
        style_row.addWidget(QLabel("Wireframe Style:"))
        self.wireframe_style_combo = QComboBox()
        self.wireframe_style_combo.addItems([
            "Single Isosurface",
            "Multiple Isosurface", 
            "Boundary Box"
        ])
        self.wireframe_style_combo.currentTextChanged.connect(self.change_wireframe_style)
        style_row.addWidget(self.wireframe_style_combo)
        wireframe_layout.addLayout(style_row)
        
        # Isosurface Level Slider with description UNDER
        self.isosurface_controls = QWidget()
        iso_layout = QVBoxLayout(self.isosurface_controls)
        
        # Slider row
        iso_slider_row = QHBoxLayout()
        iso_slider_row.addWidget(QLabel("Contour Level:"))
        self.isosurface_level_slider = QSlider(Qt.Orientation.Horizontal)
        self.isosurface_level_slider.setRange(10, 90)
        self.isosurface_level_slider.setValue(50)
        self.isosurface_level_slider.valueChanged.connect(self.update_isosurface_level)
        self.isosurface_level_label = QLabel("50% of peak")
        self.isosurface_level_label.setMinimumWidth(120)
        iso_slider_row.addWidget(self.isosurface_level_slider)
        iso_slider_row.addWidget(self.isosurface_level_label)
        iso_layout.addLayout(iso_slider_row)
        
        # Description UNDER the slider
        iso_desc = QLabel("Flux intensity level to show as wireframe contour")
        iso_desc.setStyleSheet("color: #888; font-size: 10px; margin-left: 10px;")
        iso_layout.addWidget(iso_desc)
        
        wireframe_layout.addWidget(self.isosurface_controls)
        
        # Initially hide wireframe controls
        self.wireframe_controls.setVisible(False)
        mode_layout.addRow(self.wireframe_controls)
        
        # Slice Plane Controls with description UNDER
        self.slice_controls = QWidget()
        slice_layout = QVBoxLayout(self.slice_controls)
        
        # Slice axis selector
        axis_row = QHBoxLayout()
        axis_row.addWidget(QLabel("Slice Orientation:"))
        self.slice_axis_combo = QComboBox()
        self.slice_axis_combo.addItems([
            "X-Axis (YZ Plane)", 
            "Y-Axis (XZ Plane)", 
            "Z-Axis (XY Plane)"
        ])
        self.slice_axis_combo.setCurrentText("Z-Axis (XY Plane)")
        self.slice_axis_combo.currentTextChanged.connect(self.change_slice_axis)
        axis_row.addWidget(self.slice_axis_combo)
        slice_layout.addLayout(axis_row)
        
        # Slice Position Slider with description UNDER
        slice_pos_row = QHBoxLayout()
        slice_pos_row.addWidget(QLabel("Position:"))
        self.slice_position_slider = QSlider(Qt.Orientation.Horizontal)
        self.slice_position_slider.setRange(0, 100)
        self.slice_position_slider.setValue(50)
        self.slice_position_slider.valueChanged.connect(self.update_slice_position)
        self.slice_position_label = QLabel("Center (0 km)")
        self.slice_position_label.setMinimumWidth(120)
        slice_pos_row.addWidget(self.slice_position_slider)
        slice_pos_row.addWidget(self.slice_position_label)
        slice_layout.addLayout(slice_pos_row)
        
        # Description UNDER the slider
        slice_desc = QLabel("Position of cutting plane through the 3D volume")
        slice_desc.setStyleSheet("color: #888; font-size: 10px; margin-left: 10px;")
        slice_layout.addWidget(slice_desc)
        
        # Initially hide slice controls
        self.slice_controls.setVisible(False)
        mode_layout.addRow(self.slice_controls)
        
        # Point Density Control with description UNDER
        density_row = QHBoxLayout()
        density_row.addWidget(QLabel("Point Density:"))
        self.point_density_slider = QSlider(Qt.Orientation.Horizontal)
        self.point_density_slider.setRange(500, 10000)
        self.point_density_slider.setValue(5000)
        self.point_density_slider.valueChanged.connect(self.update_point_density)
        self.point_density_label = QLabel("5,000 points")
        self.point_density_label.setMinimumWidth(100)
        density_row.addWidget(self.point_density_slider)
        density_row.addWidget(self.point_density_label)
        mode_layout.addRow(density_row)
        
        # Description UNDER
        density_desc = QLabel("Number of flux sample points to display")
        density_desc.setStyleSheet("color: #888; font-size: 10px; margin-left: 85px;")
        mode_layout.addRow("", density_desc)
        
        # Transparency control with description UNDER
        opacity_row = QHBoxLayout()
        opacity_row.addWidget(QLabel("Transparency:"))
        self.opacity_slider = QSlider(Qt.Orientation.Horizontal)
        self.opacity_slider.setRange(0, 100)
        self.opacity_slider.setValue(70)
        self.opacity_slider.valueChanged.connect(self.update_opacity)
        self.opacity_label = QLabel("70%")
        self.opacity_label.setMinimumWidth(50)
        opacity_row.addWidget(self.opacity_slider)
        opacity_row.addWidget(self.opacity_label)
        mode_layout.addRow(opacity_row)
        
        # Point Size control with description UNDER
        point_size_row = QHBoxLayout()
        point_size_row.addWidget(QLabel("Point Size:"))
        self.point_size_slider = QSlider(Qt.Orientation.Horizontal)
        self.point_size_slider.setRange(100, 800)
        self.point_size_slider.setValue(400)
        self.point_size_slider.valueChanged.connect(self.update_point_size)
        self.point_size_label = QLabel("400m radius")
        self.point_size_label.setMinimumWidth(100)
        point_size_row.addWidget(self.point_size_slider)
        point_size_row.addWidget(self.point_size_label)
        mode_layout.addRow(point_size_row)
        
        # Description UNDER
        size_desc = QLabel("Radius of each flux visualization point")
        size_desc.setStyleSheet("color: #888; font-size: 10px; margin-left: 85px;")
        mode_layout.addRow("", size_desc)
        
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
        
        viz_layout.addLayout(mode_layout)
        
        # Insert into main control layout
        if insert_index >= 0:
            control_layout.insertWidget(insert_index, viz_group)
        else:
            control_layout.addWidget(viz_group)

    def update_point_density(self, density):
        """COMPACT: Update point density with simple labeling"""
        self.point_density_label.setText(f"{density:,} points")
        
        # Only regenerate if we're in point cloud mode
        if (hasattr(self, 'viz_mode_combo') and 
            self.viz_mode_combo.currentText() == "Point Cloud" and
            hasattr(self, 'vtk_data') and self.vtk_data):
            
            if hasattr(self, '_density_update_timer'):
                self._density_update_timer.stop()
                
            self._density_update_timer = QTimer()
            self._density_update_timer.setSingleShot(True)
            self._density_update_timer.timeout.connect(lambda: self._regenerate_point_cloud(density))
            self._density_update_timer.start(200)

    def change_wireframe_style(self, style_name):
        """FIXED: Clear previous wireframe when style changes"""
        print(f"Changing wireframe style to: {style_name}")
        
        # Show/hide isosurface controls
        show_iso_controls = (style_name == "Single Isosurface")
        if hasattr(self, 'isosurface_controls'):
            self.isosurface_controls.setVisible(show_iso_controls)
        
        if self.viz_mode_combo.currentText() == "Wireframe" and self.vtk_data:
            # IMPORTANT: Clear previous wireframe first
            self.clear_field_visualization()
            
            # Create new wireframe with new style
            self.setup_wireframe_rendering()
            self.vtk_widget.GetRenderWindow().Render()

    def update_isosurface_level(self, level_percent):
        """COMPACT: Update isosurface level with concise labeling"""
        # Shorter, more concise label
        if hasattr(self, 'vtk_data') and self.vtk_data:
            scalar_range = self.vtk_data.GetScalarRange()
            actual_value = scalar_range[1] * (level_percent / 100.0)
            self.isosurface_level_label.setText(f"{level_percent}% ({actual_value:.1e})")
        else:
            self.isosurface_level_label.setText(f"{level_percent}% of peak")
        
        # Only update if we're in the right mode
        if not (self.viz_mode_combo.currentText() == "Wireframe" and 
                self.wireframe_style_combo.currentText() == "Single Isosurface" and
                self.vtk_data):
            return
        
        # FAST UPDATE: Direct contour level change
        if hasattr(self, 'field_actor') and self.field_actor:
            try:
                mapper = self.field_actor.GetMapper()
                if mapper:
                    input_conn = mapper.GetInputConnection(0, 0)
                    if input_conn:
                        edges_filter = input_conn.GetProducer()
                        contour_conn = edges_filter.GetInputConnection(0, 0)
                        if contour_conn:
                            contour_filter = contour_conn.GetProducer()
                            
                            scalar_range = self.vtk_data.GetScalarRange()
                            new_level = scalar_range[1] * (level_percent / 100.0)
                            
                            contour_filter.SetValue(0, new_level)
                            contour_filter.Modified()
                            
                            self.vtk_widget.GetRenderWindow().Render()
                            return
            except:
                pass
        
        self._schedule_isosurface_update()

    def _schedule_isosurface_update(self):
        """Schedule isosurface update with longer debounce to avoid spam"""
        if hasattr(self, '_iso_update_timer'):
            self._iso_update_timer.stop()
            
        self._iso_update_timer = QTimer()
        self._iso_update_timer.setSingleShot(True)
        self._iso_update_timer.timeout.connect(self._rebuild_isosurface)
        self._iso_update_timer.start(300)  # 300ms delay (longer to avoid spam)

    def _rebuild_isosurface(self):
        """Full rebuild of isosurface (slow method)"""
        print("Rebuilding isosurface (slow method)")
        self.clear_field_visualization()
        self.setup_wireframe_rendering()
        self.vtk_widget.GetRenderWindow().Render()
            
    def _update_isosurface(self):
        """Actually update the isosurface"""
        self.clear_field_visualization()
        self.setup_wireframe_rendering()
        self.vtk_widget.GetRenderWindow().Render()

    def change_slice_axis(self, axis_text):
        """FIXED: Clear previous slice when axis changes"""
        print(f"Changing slice axis to: {axis_text}")
        
        if self.viz_mode_combo.currentText() == "Slice Planes" and self.vtk_data:
            # IMPORTANT: Clear previous slice first
            self.clear_field_visualization()
            
            # Create new slice with new axis
            self.setup_slice_planes()
            self.vtk_widget.GetRenderWindow().Render()

    def update_slice_position(self, position_percent):
        """COMPACT: Update slice position with concise coordinate info"""
        # More compact labeling
        if hasattr(self, 'vtk_data') and self.vtk_data:
            bounds = self.vtk_data.GetBounds()
            axis_text = self.slice_axis_combo.currentText()
            
            if "X-Axis" in axis_text:
                coord_range = bounds[1] - bounds[0]
                coord_pos = bounds[0] + (position_percent/100.0) * coord_range
                axis_name = "X"
            elif "Y-Axis" in axis_text:
                coord_range = bounds[3] - bounds[2]
                coord_pos = bounds[2] + (position_percent/100.0) * coord_range
                axis_name = "Y"
            else:  # Z-Axis
                coord_range = bounds[5] - bounds[4]
                coord_pos = bounds[4] + (position_percent/100.0) * coord_range
                axis_name = "Z"
            
            # Compact label format
            if position_percent <= 25:
                position_name = "Near"
            elif position_percent >= 75:
                position_name = "Far"
            else:
                position_name = "Center"
                
            self.slice_position_label.setText(f"{position_name} ({axis_name}={coord_pos:.0f}km)")
        else:
            if position_percent <= 25:
                label = "Near"
            elif position_percent >= 75:
                label = "Far" 
            else:
                label = "Center"
            self.slice_position_label.setText(f"{label} ({position_percent}%)")
        
        # Only update if we're in slice mode
        if not (self.viz_mode_combo.currentText() == "Slice Planes" and self.vtk_data):
            return
        
        # FAST UPDATE: Direct plane origin change
        if hasattr(self, 'field_actor') and self.field_actor:
            try:
                mapper = self.field_actor.GetMapper()
                if mapper:
                    input_conn = mapper.GetInputConnection(0, 0)
                    if input_conn:
                        cutter = input_conn.GetProducer()
                        if hasattr(cutter, 'GetCutFunction'):
                            plane = cutter.GetCutFunction()
                            
                            bounds = self.vtk_data.GetBounds()
                            axis_text = self.slice_axis_combo.currentText()
                            
                            if "X-Axis" in axis_text:
                                origin_coord = bounds[0] + (position_percent/100.0) * (bounds[1] - bounds[0])
                                new_origin = [origin_coord, (bounds[2]+bounds[3])/2, (bounds[4]+bounds[5])/2]
                            elif "Y-Axis" in axis_text:
                                origin_coord = bounds[2] + (position_percent/100.0) * (bounds[3] - bounds[2])
                                new_origin = [(bounds[0]+bounds[1])/2, origin_coord, (bounds[4]+bounds[5])/2]
                            else:  # Z-Axis
                                origin_coord = bounds[4] + (position_percent/100.0) * (bounds[5] - bounds[4])
                                new_origin = [(bounds[0]+bounds[1])/2, (bounds[2]+bounds[3])/2, origin_coord]
                            
                            plane.SetOrigin(new_origin)
                            plane.Modified()
                            cutter.Modified()
                            
                            self.vtk_widget.GetRenderWindow().Render()
                            return
            except:
                pass
        
        self._schedule_slice_update()

    def _schedule_slice_update(self):
        """Schedule slice update with debounce"""
        if hasattr(self, '_slice_update_timer'):
            self._slice_update_timer.stop()
            
        self._slice_update_timer = QTimer()
        self._slice_update_timer.setSingleShot(True)
        self._slice_update_timer.timeout.connect(self._rebuild_slice)
        self._slice_update_timer.start(200)  # 200ms delay

    def _rebuild_slice(self):
        """Full rebuild of slice (slow method)"""
        print("Rebuilding slice (slow method)")
        self.clear_field_visualization()
        self.setup_slice_planes()
        self.vtk_widget.GetRenderWindow().Render()

    def _update_slice(self):
        """Actually update the slice"""
        self.clear_field_visualization()
        self.setup_slice_planes()
        self.vtk_widget.GetRenderWindow().Render()
            
    def update_point_size(self, value):
        """COMPACT: Update point size with simple labeling"""
        radius_m = value
        self.point_size_label.setText(f"{radius_m}m radius")
        
        # Only update if we're in point cloud mode and have a field actor
        if (self.viz_mode_combo.currentText() == "Point Cloud" and 
            hasattr(self, 'field_actor') and self.field_actor):
            
            self.update_point_cloud_size(radius_m)

    def _regenerate_point_cloud(self, density):
        """Regenerate point cloud with new density"""
        try:
            print(f"Regenerating point cloud with {density} points...")
            
            # Store the new target density
            self.target_point_count = density
            
            # Clear existing field visualization
            self.clear_field_visualization()
            
            # Regenerate point cloud
            self.setup_point_cloud_rendering()
            
            # Force render
            self.vtk_widget.GetRenderWindow().Render()
            
        except Exception as e:
            print(f"Error regenerating point cloud: {e}")

    def setup_point_cloud_rendering(self):
        """Enhanced point cloud with adjustable density"""
        if not self.vtk_data:
            return
            
        print("Setting up point cloud rendering...")
        
        try:
            # Get scalar data
            scalar_array = self.vtk_data.GetPointData().GetScalars()
            if not scalar_array:
                print("No scalar data for point cloud")
                return
                
            scalar_range = scalar_array.GetRange()
            num_points = self.vtk_data.GetNumberOfPoints()
            print(f"Point cloud: {num_points} points, range: {scalar_range}")
            
            # Store values
            self.current_scalar_range = scalar_range
            
            # Get target density from slider
            target_density = getattr(self, 'target_point_count', 
                                   self.point_density_slider.value() if hasattr(self, 'point_density_slider') else 5000)
            
            print(f"Target density: {target_density} points")
            
            # Analyze data distribution for smart thresholding
            print("Analyzing data distribution...")
            non_zero_count = 0
            sample_size = min(1000, scalar_array.GetSize())
            
            for i in range(sample_size):
                if scalar_array.GetValue(i) > 0:
                    non_zero_count += 1
                    
            non_zero_fraction = non_zero_count / sample_size
            print(f"Non-zero values: {non_zero_count}/{sample_size} ({non_zero_fraction:.1%})")
            
            # Adaptive threshold
            if non_zero_fraction > 0.5:
                threshold_fraction = 0.02  # 2% for dense data
            elif non_zero_fraction > 0.1:
                threshold_fraction = 0.005  # 0.5% for medium data
            else:
                threshold_fraction = 0.001  # 0.1% for sparse data
            
            # Apply threshold
            threshold = vtk.vtkThreshold()
            threshold.SetInputData(self.vtk_data)
            threshold_value = scalar_range[1] * threshold_fraction
            threshold.SetLowerThreshold(threshold_value)
            threshold.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_UPPER)
            threshold.Update()
            
            significant_data = threshold.GetOutput()
            significant_points = significant_data.GetNumberOfPoints()
            print(f"Significant points (>{threshold_value:.2e}): {significant_points}")
            
            # Fallback to lower thresholds if needed
            if significant_points == 0:
                for lower_fraction in [0.001, 0.0001, 0.0]:
                    threshold_value = scalar_range[1] * lower_fraction if lower_fraction > 0 else 0
                    threshold.SetLowerThreshold(threshold_value)
                    threshold.Update()
                    significant_data = threshold.GetOutput()
                    significant_points = significant_data.GetNumberOfPoints()
                    if significant_points > 0 or lower_fraction == 0:
                        break
                        
                if significant_points == 0:
                    significant_data = self.vtk_data
                    significant_points = num_points
            
            # Subsample to target density
            if significant_points > target_density:
                print(f"Subsampling {significant_points} -> {target_density} points")
                mask = vtk.vtkMaskPoints()
                mask.SetInputData(significant_data)
                mask.SetMaximumNumberOfPoints(target_density)
                mask.SetRandomMode(True)
                mask.Update()
                final_data = mask.GetOutput()
            else:
                final_data = significant_data
                
            final_count = final_data.GetNumberOfPoints()
            print(f"Final point count: {final_count}")
            
            if final_count == 0:
                print("ERROR: No points to visualize!")
                return
            
            # Add jitter and create visualization
            jittered_data = self.add_spatial_jitter(final_data)
            self.current_final_data = jittered_data
            
            # Create glyphs
            initial_radius = self.point_size_slider.value() if hasattr(self, 'point_size_slider') else 400
            self.create_point_cloud_glyphs(jittered_data, initial_radius)
            
            print(f"Point cloud complete: {final_count} spheres")
            
        except Exception as e:
            print(f"Point cloud rendering failed: {e}")
            import traceback
            traceback.print_exc()

    def update_point_size(self, value):
        """Update point size for point cloud visualization"""
        radius_m = value  # Slider value is in meters
        self.point_size_label.setText(f"{radius_m}m")
        
        # Only update if we're in point cloud mode and have a field actor
        if (self.viz_mode_combo.currentText() == "Point Cloud" and 
            hasattr(self, 'field_actor') and self.field_actor):
            
            print(f"Updating point size to {radius_m}m radius")
            
            # Re-create the point cloud with new size
            self.update_point_cloud_size(radius_m)

    def update_point_cloud_size(self, radius_m):
        """Fast point cloud size update - OPTIMIZED"""
        if not hasattr(self, 'current_final_data') or not self.current_final_data:
            return
            
        # Use a timer to debounce rapid slider changes
        if hasattr(self, '_size_update_timer'):
            self._size_update_timer.stop()
            
        self._size_update_timer = QTimer()
        self._size_update_timer.setSingleShot(True)
        self._size_update_timer.timeout.connect(lambda: self._do_size_update(radius_m))
        self._size_update_timer.start(50)  # 50ms delay to debounce

    def _do_size_update(self, radius_m):
        """Actually perform the size update"""
        try:
            print(f"Updating to radius {radius_m}m...")
            
            # Quick method: just update the sphere source if possible
            if (hasattr(self, 'field_actor') and self.field_actor and 
                hasattr(self, 'current_sphere_source')):
                
                # Try to update existing sphere source
                self.current_sphere_source.SetRadius(radius_m)
                self.current_sphere_source.Modified()
                
            else:
                # Fallback: recreate glyphs (slower but more reliable)
                self.create_point_cloud_glyphs(self.current_final_data, radius_m)
            
            # Force render
            self.vtk_widget.GetRenderWindow().Render()
            
        except Exception as e:
            print(f"Error in size update: {e}")
            # Fallback to full recreation
            self.create_point_cloud_glyphs(self.current_final_data, radius_m)
            self.vtk_widget.GetRenderWindow().Render()

    def update_threshold_display(self):
        """Update threshold labels without applying"""
        min_val = self.threshold_min_slider.value()
        max_val = self.threshold_max_slider.value()
        self.threshold_min_label.setText(f"{min_val}%")
        self.threshold_max_label.setText(f"{max_val}%")

    def apply_threshold_filter(self):
        """Apply threshold filter to current visualization"""
        if self.viz_mode_combo.currentText() == "Point Cloud":
            # For point cloud, just update with current threshold
            current_size = self.point_size_slider.value()
            self.update_point_cloud_size(current_size)
        else:
            # For other modes, use existing threshold logic
            self.update_threshold()
            
    def change_visualization_mode(self, mode):
        """Enhanced mode change with control visibility"""
        # Show/hide relevant controls based on mode
        self.wireframe_controls.setVisible(mode == "Wireframe")
        self.slice_controls.setVisible(mode == "Slice Planes")
        
        # Call the original visualization change
        if not self.vtk_data:
            return
        
        print(f"=== CHANGING VISUALIZATION MODE TO: {mode} ===")
        
        try:
            self.clear_field_visualization()
            
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
                self.setup_field_visualization()
                
            self.vtk_widget.GetRenderWindow().Render()
            print(f"=== VISUALIZATION MODE CHANGED TO: {mode} ===")
        
        except Exception as e:
            print(f"Error changing visualization mode: {e}")
            self.setup_field_visualization()
            self.vtk_widget.GetRenderWindow().Render()

    def clear_field_visualization(self):
        """ENHANCED clearing to handle ALL visualization actors properly"""
        print("=== CLEARING ALL FIELD VISUALIZATION ===")
        
        actors_removed = 0
        
        # Remove field actor (used by most modes)
        if hasattr(self, 'field_actor') and self.field_actor:
            self.renderer.RemoveActor(self.field_actor)
            self.field_actor = None
            actors_removed += 1
            print("Removed field_actor")
            
        # Remove volume actor
        if hasattr(self, 'volume_actor') and self.volume_actor:
            self.renderer.RemoveVolume(self.volume_actor)
            self.volume_actor = None
            actors_removed += 1
            print("Removed volume_actor")
            
        # Remove wireframe actors (for multiple wireframe mode)
        if hasattr(self, 'wireframe_actors'):
            for i, actor in enumerate(self.wireframe_actors):
                if actor:
                    self.renderer.RemoveActor(actor)
                    actors_removed += 1
            self.wireframe_actors = []
            print(f"Removed wireframe_actors")
            
        # Remove slice actors (for multiple slice mode)
        if hasattr(self, 'slice_actors'):
            for i, actor in enumerate(self.slice_actors):
                if actor:
                    self.renderer.RemoveActor(actor)
                    actors_removed += 1
            self.slice_actors = []
            print(f"Removed slice_actors")
            
        # Remove isosurface actors
        if hasattr(self, 'isosurface_actors'):
            for i, actor in enumerate(self.isosurface_actors):
                if actor:
                    self.renderer.RemoveActor(actor)
                    actors_removed += 1
            self.isosurface_actors = []
            print(f"Removed isosurface_actors")
            
        # IMPORTANT: Remove scalar bar (it persists across modes)
        if hasattr(self, 'scalar_bar') and self.scalar_bar:
            self.renderer.RemoveViewProp(self.scalar_bar)
            self.scalar_bar = None
            print("Removed scalar_bar")
            
        print(f"Total actors removed: {actors_removed}")
        print("=== CLEAR COMPLETE ===")

    def setup_volume_rendering(self):
        """Fixed volume rendering with better error handling"""
        if not self.vtk_data:
            return
            
        print("Setting up volume rendering...")
        
        try:
            # Check if we already have image data
            if isinstance(self.vtk_data, vtk.vtkImageData):
                volume_data = self.vtk_data
                print("Using existing image data")
            else:
                # Convert to image data using resample filter
                print("Converting to image data for volume rendering...")
                
                # Get bounds
                bounds = self.vtk_data.GetBounds()
                
                # Create a reasonable grid resolution
                spacing = [(bounds[1]-bounds[0])/30, (bounds[3]-bounds[2])/30, (bounds[5]-bounds[4])/30]
                
                # Resample to structured grid
                resample = vtk.vtkResampleToImage()
                resample.SetInputData(self.vtk_data)
                resample.SetSamplingDimensions(30, 30, 30)  # Reduced for performance
                resample.Update()
                volume_data = resample.GetOutput()
                
                if volume_data.GetNumberOfPoints() == 0:
                    print("Failed to create volume data")
                    return
                    
            # Use CPU volume mapper for stability
            volume_mapper = vtk.vtkFixedPointVolumeRayCastMapper()
            volume_mapper.SetInputData(volume_data)
            volume_mapper.SetBlendModeToComposite()
            
            # Volume properties
            volume_property = vtk.vtkVolumeProperty()
            volume_property.SetInterpolationTypeToLinear()
            volume_property.ShadeOn()
            
            # Get scalar range
            scalar_range = volume_data.GetScalarRange()
            print(f"Volume scalar range: {scalar_range}")
            
            if scalar_range[1] > scalar_range[0]:
                # Color transfer function
                color_func = vtk.vtkColorTransferFunction()
                color_func.AddRGBPoint(scalar_range[0], 0.0, 0.0, 0.2)
                color_func.AddRGBPoint(scalar_range[1] * 0.25, 0.0, 0.0, 1.0)
                color_func.AddRGBPoint(scalar_range[1] * 0.5, 0.0, 1.0, 0.0)
                color_func.AddRGBPoint(scalar_range[1] * 0.75, 1.0, 1.0, 0.0)
                color_func.AddRGBPoint(scalar_range[1], 1.0, 0.0, 0.0)
                
                # Opacity transfer function - make more transparent
                opacity_func = vtk.vtkPiecewiseFunction()
                opacity_func.AddPoint(scalar_range[0], 0.0)
                opacity_func.AddPoint(scalar_range[1] * 0.1, 0.0)
                opacity_func.AddPoint(scalar_range[1] * 0.3, 0.05)
                opacity_func.AddPoint(scalar_range[1] * 0.6, 0.1)
                opacity_func.AddPoint(scalar_range[1], 0.2)
                
                volume_property.SetColor(color_func)
                volume_property.SetScalarOpacity(opacity_func)
                
            # Create volume
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
        """Fixed isosurface rendering"""
        if not self.vtk_data:
            return
            
        print("Setting up isosurface rendering...")
        
        try:
            # Get scalar range
            scalar_range = self.vtk_data.GetScalarRange()
            print(f"Scalar range for isosurfaces: {scalar_range}")
            
            if scalar_range[1] <= scalar_range[0]:
                print("Invalid scalar range for isosurfaces")
                return
                
            # Create contour filter
            contour = vtk.vtkContourFilter()
            contour.SetInputData(self.vtk_data)
            
            # Generate isosurface values - only in meaningful range
            num_contours = self.num_isosurfaces_spin.value()
            min_val = scalar_range[1] * 0.3  # Start at 30% of max
            max_val = scalar_range[1] * 0.9  # End at 90% of max
            
            if max_val > min_val:
                contour.GenerateValues(num_contours, min_val, max_val)
                contour.Update()
                
                output = contour.GetOutput()
                print(f"Contour generated {output.GetNumberOfPoints()} points")
                
                if output.GetNumberOfPoints() > 0:
                    # Create mapper and actor
                    contour_mapper = vtk.vtkPolyDataMapper()
                    contour_mapper.SetInputConnection(contour.GetOutputPort())
                    contour_mapper.SetScalarRange(scalar_range)
                    
                    # Setup lookup table
                    lut = self.create_lookup_table(self.colormap_combo.currentText())
                    contour_mapper.SetLookupTable(lut)
                    
                    self.field_actor = vtk.vtkActor()
                    self.field_actor.SetMapper(contour_mapper)
                    self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
                    
                    self.renderer.AddActor(self.field_actor)
                    print(f"Created {num_contours} isosurfaces")
                else:
                    print("No isosurfaces generated - trying point cloud fallback")
                    self.setup_point_cloud_rendering()
            else:
                print("Invalid contour range")
                
        except Exception as e:
            print(f"Isosurface rendering failed: {e}")
            self.setup_point_cloud_rendering()

    # def setup_point_cloud_rendering(self):
    #     """Fixed point cloud with better threshold logic"""
    #     if not self.vtk_data:
    #         return
            
    #     print("Setting up optimized point cloud rendering...")
        
    #     try:
    #         # Get scalar data
    #         scalar_array = self.vtk_data.GetPointData().GetScalars()
    #         if not scalar_array:
    #             print("No scalar data for point cloud")
    #             return
                
    #         scalar_range = scalar_array.GetRange()
    #         num_points = self.vtk_data.GetNumberOfPoints()
    #         print(f"Point cloud: {num_points} points, range: {scalar_range}")
            
    #         # Store these values for later use
    #         self.current_scalar_range = scalar_range
            
    #         # Better threshold logic - check the data distribution first
    #         print("Analyzing data distribution...")
    #         non_zero_count = 0
    #         sample_size = min(1000, scalar_array.GetSize())
            
    #         for i in range(sample_size):
    #             if scalar_array.GetValue(i) > 0:
    #                 non_zero_count += 1
                    
    #         non_zero_fraction = non_zero_count / sample_size
    #         print(f"Non-zero values in sample: {non_zero_count}/{sample_size} ({non_zero_fraction:.1%})")
            
    #         # Adaptive threshold based on data distribution
    #         if non_zero_fraction > 0.5:
    #             # Most data is non-zero, use higher threshold
    #             threshold_fraction = 0.05  # 5%
    #         elif non_zero_fraction > 0.1:
    #             # Some data is non-zero, use medium threshold  
    #             threshold_fraction = 0.01  # 1%
    #         else:
    #             # Most data is zero, use very low threshold
    #             threshold_fraction = 0.001  # 0.1%
            
    #         self.original_threshold_fraction = threshold_fraction
            
    #         # Apply threshold
    #         threshold = vtk.vtkThreshold()
    #         threshold.SetInputData(self.vtk_data)
    #         threshold_value = scalar_range[1] * threshold_fraction
    #         print(f"Using threshold: {threshold_value:.2e} ({threshold_fraction:.1%} of max)")
            
    #         threshold.SetLowerThreshold(threshold_value)
    #         threshold.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_UPPER)
    #         threshold.Update()
            
    #         significant_data = threshold.GetOutput()
    #         significant_points = significant_data.GetNumberOfPoints()
    #         print(f"Significant flux points: {significant_points}")
            
    #         # If still no points, try progressively lower thresholds
    #         if significant_points == 0:
    #             print("No points found, trying lower thresholds...")
                
    #             for lower_fraction in [0.001, 0.0001, 0.00001, 0.0]:
    #                 threshold_value = scalar_range[1] * lower_fraction
    #                 threshold.SetLowerThreshold(threshold_value)
    #                 threshold.Update()
    #                 significant_data = threshold.GetOutput()
    #                 significant_points = significant_data.GetNumberOfPoints()
    #                 print(f"  Threshold {threshold_value:.2e}: {significant_points} points")
                    
    #                 if significant_points > 0:
    #                     self.original_threshold_fraction = lower_fraction
    #                     break
                        
    #             # If still nothing, use all data
    #             if significant_points == 0:
    #                 print("Using all data (no threshold)")
    #                 significant_data = self.vtk_data
    #                 significant_points = num_points
    #                 self.original_threshold_fraction = 0.0
            
    #         # Subsample for performance
    #         target_points = 2000
    #         if significant_points > target_points:
    #             print(f"Subsampling to {target_points} points for performance...")
    #             mask = vtk.vtkMaskPoints()
    #             mask.SetInputData(significant_data)
    #             mask.SetMaximumNumberOfPoints(target_points)
    #             mask.SetRandomMode(True)
    #             mask.Update()
    #             final_data = mask.GetOutput()
    #         else:
    #             final_data = significant_data
                
    #         final_point_count = final_data.GetNumberOfPoints()
    #         print(f"Final point count: {final_point_count}")
            
    #         if final_point_count == 0:
    #             print("ERROR: No points to visualize!")
    #             return
            
    #         # Add jitter to break up grid artifacts
    #         jittered_data = self.add_spatial_jitter(final_data)
            
    #         # Store for later updates
    #         self.current_final_data = jittered_data
            
    #         # Create initial visualization
    #         initial_radius = self.point_size_slider.value() if hasattr(self, 'point_size_slider') else 400
    #         self.create_point_cloud_glyphs(jittered_data, initial_radius)
            
    #         print(f"Point cloud setup complete: {jittered_data.GetNumberOfPoints()} visible spheres")
            
    #     except Exception as e:
    #         print(f"Point cloud rendering failed: {e}")
    #         import traceback
    #         traceback.print_exc()

    def add_spatial_jitter(self, input_data):
        """Add random spatial jitter to break up grid artifacts"""
        if input_data.GetNumberOfPoints() == 0:
            print("No points to jitter")
            return input_data
            
        print("Adding spatial jitter to reduce grid artifacts...")
        
        try:
            # Get the spacing of the original structured grid
            jitter_radius = 300.0  # km
            
            # Create a new points array with jitter
            original_points = input_data.GetPoints()
            n_points = original_points.GetNumberOfPoints()
            
            jittered_points = vtk.vtkPoints()
            
            import random
            random.seed(42)  # Reproducible jitter
            
            for i in range(n_points):
                # Get original point
                orig_point = original_points.GetPoint(i)
                
                # Add random jitter in all three dimensions
                jitter_x = (random.random() - 0.5) * 2 * jitter_radius
                jitter_y = (random.random() - 0.5) * 2 * jitter_radius  
                jitter_z = (random.random() - 0.5) * 2 * jitter_radius
                
                # Create jittered point
                new_point = [
                    orig_point[0] + jitter_x,
                    orig_point[1] + jitter_y,
                    orig_point[2] + jitter_z
                ]
                
                jittered_points.InsertNextPoint(new_point)
            
            # Create new polydata with jittered points
            jittered_polydata = vtk.vtkPolyData()
            jittered_polydata.SetPoints(jittered_points)
            
            # Copy scalar data
            jittered_polydata.GetPointData().SetScalars(input_data.GetPointData().GetScalars())
            
            # Create vertices
            vertices = vtk.vtkCellArray()
            for i in range(n_points):
                vertex = vtk.vtkVertex()
                vertex.GetPointIds().SetId(0, i)
                vertices.InsertNextCell(vertex)
            jittered_polydata.SetVerts(vertices)
            
            print(f"Applied jitter to {n_points} points (±{jitter_radius}km)")
            return jittered_polydata
            
        except Exception as e:
            print(f"Error adding jitter: {e}")
            return input_data  # Return original if jitter fails

    def create_point_cloud_glyphs(self, point_data, radius):
        """Create the actual glyph visualization"""
        if point_data.GetNumberOfPoints() == 0:
            print("No points to create glyphs for")
            return
            
        try:
            print(f"Creating glyphs for {point_data.GetNumberOfPoints()} points...")
            
            # Use simpler spheres for better performance
            sphere = vtk.vtkSphereSource()
            sphere.SetRadius(radius)
            sphere.SetThetaResolution(8)  # Reduced for performance
            sphere.SetPhiResolution(8)    # Reduced for performance
            
            # Store sphere source for potential fast updates
            self.current_sphere_source = sphere
            
            # Create glyph
            glyph = vtk.vtkGlyph3D()
            glyph.SetInputData(point_data)
            glyph.SetSourceConnection(sphere.GetOutputPort())
            glyph.SetScaleModeToDataScalingOff()
            glyph.SetColorModeToColorByScalar()
            glyph.Update()
            
            print(f"Glyph filter created {glyph.GetOutput().GetNumberOfPoints()} glyph points")
            
            # Create mapper
            glyph_mapper = vtk.vtkPolyDataMapper()
            glyph_mapper.SetInputConnection(glyph.GetOutputPort())
            glyph_mapper.SetScalarRange(self.current_scalar_range)
            glyph_mapper.ScalarVisibilityOn()
            
            # Setup lookup table
            lut = self.create_lookup_table(self.colormap_combo.currentText() if hasattr(self, 'colormap_combo') else 'Blue to Red')
            glyph_mapper.SetLookupTable(lut)
            
            # Create or update actor
            if hasattr(self, 'field_actor') and self.field_actor:
                # Update existing actor
                self.field_actor.SetMapper(glyph_mapper)
            else:
                # Create new actor
                self.field_actor = vtk.vtkActor()
                self.field_actor.SetMapper(glyph_mapper)
                self.renderer.AddActor(self.field_actor)
            
            # Set opacity
            opacity = self.opacity_slider.value() / 100.0 if hasattr(self, 'opacity_slider') else 0.7
            self.field_actor.GetProperty().SetOpacity(opacity)
            
            # Setup scalar bar
            if hasattr(self, 'current_scalar_range'):
                scalar_name = self.vtk_data.GetPointData().GetScalars().GetName()
                self.setup_scalar_bar(lut, scalar_name)
            
            print(f"Glyphs created successfully: {point_data.GetNumberOfPoints()} spheres, radius {radius}m")
            
        except Exception as e:
            print(f"Error creating glyphs: {e}")
            import traceback
            traceback.print_exc()

    def setup_simple_points(self):
        """Fallback: simple point rendering without glyphs"""
        print("Setting up simple point rendering...")
        
        try:
            # Simple mapper without glyphs
            point_mapper = vtk.vtkDataSetMapper()
            point_mapper.SetInputData(self.vtk_data)
            
            scalar_range = self.vtk_data.GetScalarRange()
            point_mapper.SetScalarRange(scalar_range)
            
            lut = self.create_lookup_table(self.colormap_combo.currentText())
            point_mapper.SetLookupTable(lut)
            
            self.field_actor = vtk.vtkActor()
            self.field_actor.SetMapper(point_mapper)
            self.field_actor.GetProperty().SetPointSize(4.0)
            self.field_actor.GetProperty().SetRenderPointsAsSpheres(True)
            self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
            
            self.renderer.AddActor(self.field_actor)
            print("Simple point rendering setup complete")
            
        except Exception as e:
            print(f"Simple point rendering failed: {e}")

    def setup_wireframe_rendering(self):
        """Enhanced wireframe with user-controllable options"""
        if not self.vtk_data:
            return
        
        print("Setting up enhanced wireframe rendering...")
        
        try:
            scalar_array = self.vtk_data.GetPointData().GetScalars()
            if not scalar_array:
                print("No scalar data for wireframe")
                return
            
            scalar_range = scalar_array.GetRange()
        
            # Get current wireframe style
            style = getattr(self, 'wireframe_style_combo', None)
            if style:
                current_style = style.currentText()
            else:
                current_style = "Single Isosurface"
                
            if current_style == "Single Isosurface":
                self.create_controllable_single_isosurface(scalar_range)
            elif current_style == "Multiple Isosurface":
                self.create_multiple_isosurface_wireframes(scalar_range)
            else:  # Boundary Box
                self.create_boundary_wireframe()
            
            print("Enhanced wireframe rendering complete")
        
        except Exception as e:
            print(f"Wireframe rendering failed: {e}")

    def create_controllable_single_isosurface(self, scalar_range):
        """Create isosurface with stored references for fast updates"""
        if scalar_range[1] <= scalar_range[0]:
            return
            
        level_percent = getattr(self, 'isosurface_level_slider', None)
        if level_percent:
            percent = level_percent.value() / 100.0
        else:
            percent = 0.5
            
        contour_level = scalar_range[1] * percent
        print(f"Creating isosurface wireframe at {percent*100:.0f}% ({contour_level:.2e})")
        
        # Create contour filter
        self.current_contour_filter = vtk.vtkContourFilter()  # Store reference for fast updates
        self.current_contour_filter.SetInputData(self.vtk_data)
        self.current_contour_filter.SetValue(0, contour_level)
        self.current_contour_filter.Update()
        
        contour_output = self.current_contour_filter.GetOutput()
        
        if contour_output.GetNumberOfPoints() > 0:
            # Create edges filter
            self.current_edges_filter = vtk.vtkExtractEdges()  # Store reference
            self.current_edges_filter.SetInputConnection(self.current_contour_filter.GetOutputPort())
            self.current_edges_filter.Update()
            
            wireframe_mapper = vtk.vtkPolyDataMapper()
            wireframe_mapper.SetInputConnection(self.current_edges_filter.GetOutputPort())
            
            self.field_actor = vtk.vtkActor()
            self.field_actor.SetMapper(wireframe_mapper)
            self.field_actor.GetProperty().SetColor(0.9, 0.9, 0.2)
            self.field_actor.GetProperty().SetLineWidth(2.0)
            self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
            
            self.renderer.AddActor(self.field_actor)
            print(f"Isosurface wireframe created: {contour_output.GetNumberOfPoints()} points")
        else:
            print("No contour generated at this level")

    def create_isosurface_wireframes(self, scalar_range):
        """Create wireframe contours at different flux levels"""
        if scalar_range[1] <= scalar_range[0]:
            return
        
        print("Creating isosurface wireframes...")
    
        # Create 3-4 isosurface levels in the meaningful range
        min_val = scalar_range[1] * 0.2   # 20% of max
        max_val = scalar_range[1] * 0.8   # 80% of max
    
        if max_val > min_val:
            num_levels = 4
            levels = np.linspace(min_val, max_val, num_levels)
        
            for i, level in enumerate(levels):
                # Create contour
                contour = vtk.vtkContourFilter()
                contour.SetInputData(self.vtk_data)
                contour.SetValue(0, level)
                contour.Update()
            
                contour_output = contour.GetOutput()
            
                if contour_output.GetNumberOfPoints() > 0:
                    # Extract edges to create wireframe
                    edges = vtk.vtkExtractEdges()
                    edges.SetInputData(contour_output)
                    edges.Update()
                
                    # Create mapper
                    wireframe_mapper = vtk.vtkPolyDataMapper()
                    wireframe_mapper.SetInputConnection(edges.GetOutputPort())
                
                    # Create actor
                    wireframe_actor = vtk.vtkActor()
                    wireframe_actor.SetMapper(wireframe_mapper)
                
                    # Color by level (blue to red)
                    intensity = i / (num_levels - 1)
                    color = [intensity, 0.5, 1.0 - intensity]
                    wireframe_actor.GetProperty().SetColor(color)
                    wireframe_actor.GetProperty().SetLineWidth(2.0)
                    wireframe_actor.GetProperty().SetOpacity(0.8)
                
                    self.wireframe_actors.append(wireframe_actor)
                    self.renderer.AddActor(wireframe_actor)
                
                    print(f"  Level {i+1}: {level:.2e} ({contour_output.GetNumberOfPoints()} points)")

    def create_slice_wireframes(self):
        """Create wireframe on cross-sectional planes"""
        if not self.vtk_data:
            return
        
        print("Adding slice wireframes...")
    
        try:
            # Get bounds for slice positioning
            bounds = self.vtk_data.GetBounds()
            center = [(bounds[1]+bounds[0])/2, (bounds[3]+bounds[2])/2, (bounds[5]+bounds[4])/2]
        
            # Create one slice plane through the center
            plane = vtk.vtkPlane()
            plane.SetOrigin(center)
            plane.SetNormal(0, 0, 1)  # XY plane
        
            # Cut the data
            cutter = vtk.vtkCutter()
            cutter.SetInputData(self.vtk_data)
            cutter.SetCutFunction(plane)
            cutter.Update()
        
            slice_data = cutter.GetOutput()
        
            if slice_data.GetNumberOfPoints() > 0:
                # Extract edges
                edges = vtk.vtkExtractEdges()
                edges.SetInputData(slice_data)
                edges.Update()
            
                # Create mapper
                slice_wireframe_mapper = vtk.vtkPolyDataMapper()
                slice_wireframe_mapper.SetInputConnection(edges.GetOutputPort())
                
                # Setup color mapping
                scalar_range = self.vtk_data.GetScalarRange()
                slice_wireframe_mapper.SetScalarRange(scalar_range)
            
                lut = self.create_lookup_table(self.colormap_combo.currentText())
                slice_wireframe_mapper.SetLookupTable(lut)
            
                # Create actor
                slice_wireframe_actor = vtk.vtkActor()
                slice_wireframe_actor.SetMapper(slice_wireframe_mapper)
                slice_wireframe_actor.GetProperty().SetLineWidth(1.5)
                slice_wireframe_actor.GetProperty().SetOpacity(0.6)
            
                self.wireframe_actors.append(slice_wireframe_actor)
                self.renderer.AddActor(slice_wireframe_actor)
            
                print(f"  Added slice wireframe ({slice_data.GetNumberOfPoints()} points)")
            
        except Exception as e:
            print(f"Error creating slice wireframes: {e}")

    def create_simple_boundary_wireframe(self):
        """Fallback: create simple boundary wireframe (your original approach)"""
        print("Creating simple boundary wireframe...")
    
        try:
            # This is your original approach - kept as fallback
            if isinstance(self.vtk_data, vtk.vtkUnstructuredGrid):
                edges = vtk.vtkExtractEdges()
                edges.SetInputData(self.vtk_data)
                edges.Update()
                wireframe_data = edges.GetOutput()
            else:
                # For other types, extract surface then edges
                surface = vtk.vtkDataSetSurfaceFilter()
                surface.SetInputData(self.vtk_data)
                surface.Update()
                
                edges = vtk.vtkExtractEdges()
                edges.SetInputData(surface.GetOutput())
                edges.Update()
                wireframe_data = edges.GetOutput()
                
            print(f"Simple wireframe: {wireframe_data.GetNumberOfCells()} edges")
        
            # Create mapper
            wireframe_mapper = vtk.vtkPolyDataMapper()
            wireframe_mapper.SetInputData(wireframe_data)
            
            scalar_range = self.vtk_data.GetScalarRange()
            wireframe_mapper.SetScalarRange(scalar_range)
            
            lut = self.create_lookup_table(self.colormap_combo.currentText())
            wireframe_mapper.SetLookupTable(lut)
            
            self.field_actor = vtk.vtkActor()
            self.field_actor.SetMapper(wireframe_mapper)
            self.field_actor.GetProperty().SetLineWidth(1.0)
            self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
            
            self.renderer.AddActor(self.field_actor)
            print("Simple boundary wireframe created")
        
        except Exception as e:
            print(f"Simple boundary wireframe failed: {e}")

    def setup_surface_with_edges(self):
        """Fixed surface with edges rendering"""
        if not self.vtk_data:
            return
            
        print("Setting up surface with edges...")
        
        try:
            # Extract outer surface
            surface = vtk.vtkDataSetSurfaceFilter()
            surface.SetInputData(self.vtk_data)
            surface.Update()
            
            surface_data = surface.GetOutput()
            print(f"Extracted surface with {surface_data.GetNumberOfPoints()} points")
            
            if surface_data.GetNumberOfPoints() == 0:
                print("No surface extracted - falling back to point cloud")
                self.setup_point_cloud_rendering()
                return
            
            # Create mapper
            surface_mapper = vtk.vtkPolyDataMapper()
            surface_mapper.SetInputData(surface_data)
            
            scalar_range = self.vtk_data.GetScalarRange()
            surface_mapper.SetScalarRange(scalar_range)
            
            lut = self.create_lookup_table(self.colormap_combo.currentText())
            surface_mapper.SetLookupTable(lut)
            
            self.field_actor = vtk.vtkActor()
            self.field_actor.SetMapper(surface_mapper)
            self.field_actor.GetProperty().SetRepresentationToSurface()
            self.field_actor.GetProperty().EdgeVisibilityOn()
            self.field_actor.GetProperty().SetEdgeColor(0.2, 0.2, 0.2)
            self.field_actor.GetProperty().SetLineWidth(1.0)
            self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
            
            self.renderer.AddActor(self.field_actor)
            print("Surface with edges setup complete")
            
        except Exception as e:
            print(f"Surface with edges failed: {e}")

    def setup_slice_planes(self):
        """Create slice with stored references for fast updates"""
        if not self.vtk_data:
            return
            
        print("Setting up enhanced slice planes...")
        
        try:
            bounds = self.vtk_data.GetBounds()
            
            # Get current settings
            axis_combo = getattr(self, 'slice_axis_combo', None)
            axis_text = axis_combo.currentText() if axis_combo else "Z-Axis (XY Plane)"
            
            pos_slider = getattr(self, 'slice_position_slider', None)
            position_percent = pos_slider.value() / 100.0 if pos_slider else 0.5
            
            # Determine axis and position
            if "X-Axis" in axis_text:
                normal = [1, 0, 0]
                origin_coord = bounds[0] + position_percent * (bounds[1] - bounds[0])
                origin = [origin_coord, (bounds[2]+bounds[3])/2, (bounds[4]+bounds[5])/2]
                axis_name = "X"
            elif "Y-Axis" in axis_text:
                normal = [0, 1, 0]
                origin_coord = bounds[2] + position_percent * (bounds[3] - bounds[2])
                origin = [(bounds[0]+bounds[1])/2, origin_coord, (bounds[4]+bounds[5])/2]
                axis_name = "Y"
            else:  # Z-Axis
                normal = [0, 0, 1]
                origin_coord = bounds[4] + position_percent * (bounds[5] - bounds[4])
                origin = [(bounds[0]+bounds[1])/2, (bounds[2]+bounds[3])/2, origin_coord]
                axis_name = "Z"
                
            print(f"Creating {axis_name}-axis slice at {position_percent*100:.0f}% ({origin_coord:.0f} km)")
            
            # Create cutting plane (store reference for fast updates)
            self.current_slice_plane = vtk.vtkPlane()
            self.current_slice_plane.SetOrigin(origin)
            self.current_slice_plane.SetNormal(normal)
            
            # Create cutter (store reference for fast updates)
            self.current_slice_cutter = vtk.vtkCutter()
            self.current_slice_cutter.SetInputData(self.vtk_data)
            self.current_slice_cutter.SetCutFunction(self.current_slice_plane)
            self.current_slice_cutter.Update()
            
            slice_data = self.current_slice_cutter.GetOutput()
            
            if slice_data.GetNumberOfPoints() > 0:
                slice_mapper = vtk.vtkPolyDataMapper()
                slice_mapper.SetInputConnection(self.current_slice_cutter.GetOutputPort())
                
                scalar_range = self.vtk_data.GetScalarRange()
                slice_mapper.SetScalarRange(scalar_range)
                
                lut = self.create_lookup_table(self.colormap_combo.currentText())
                slice_mapper.SetLookupTable(lut)
                
                self.field_actor = vtk.vtkActor()
                self.field_actor.SetMapper(slice_mapper)
                self.field_actor.GetProperty().SetOpacity(self.opacity_slider.value() / 100.0)
                
                self.renderer.AddActor(self.field_actor)
                
                # Setup scalar bar
                self.setup_scalar_bar(lut, "electron_flux")
                
                print(f"Solid slice plane created: {slice_data.GetNumberOfPoints()} points")
            else:
                print("No slice data generated")
                
        except Exception as e:
            print(f"Slice planes failed: {e}")

    def create_lookup_table(self, colormap_name):
        """Create lookup table with better color schemes"""
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(256)
        
        try:
            if colormap_name == "Blue to Red":
                lut.SetHueRange(0.667, 0.0)  # Blue to red
                lut.SetSaturationRange(1.0, 1.0)
                lut.SetValueRange(0.3, 1.0)
            elif colormap_name == "Viridis":
                # Approximate viridis colormap
                for i in range(256):
                    t = i / 255.0
                    if t < 0.25:
                        r, g, b = 0.267004 + t*0.183371, 0.004874 + t*0.478894, 0.329415 + t*0.511095
                    elif t < 0.5:
                        t_norm = (t - 0.25) / 0.25
                        r, g, b = 0.127568 + t_norm*0.253935, 0.566949 + t_norm*0.214982, 0.550556 + t_norm*(-0.133721)
                    elif t < 0.75:
                        t_norm = (t - 0.5) / 0.25
                        r, g, b = 0.369214 + t_norm*0.304662, 0.788675 + t_norm*0.138379, 0.382914 + t_norm*(-0.224125)
                    else:
                        t_norm = (t - 0.75) / 0.25
                        r, g, b = 0.663765 + t_norm*0.267004, 0.865006 + t_norm*0.104775, 0.197275 + t_norm*0.109709
                    lut.SetTableValue(i, r, g, b, 1.0)
            elif colormap_name == "Plasma":
                # Approximate plasma colormap
                for i in range(256):
                    t = i / 255.0
                    r = 0.050383 + t * (0.940015 - 0.050383)
                    g = 0.029803 + t * t * (0.975158 - 0.029803)
                    b = 0.527975 + t * (0.131326 - 0.527975)
                    lut.SetTableValue(i, r, g, b, 1.0)
            elif colormap_name == "Cool to Warm":
                lut.SetHueRange(0.667, 0.0)
                lut.SetSaturationRange(1.0, 1.0)
                lut.SetValueRange(0.4, 1.0)
            elif colormap_name == "Rainbow":
                lut.SetHueRange(0.0, 0.667)
                lut.SetSaturationRange(1.0, 1.0)
                lut.SetValueRange(1.0, 1.0)
            elif colormap_name == "Grayscale":
                lut.SetHueRange(0.0, 0.0)
                lut.SetSaturationRange(0.0, 0.0)
                lut.SetValueRange(0.0, 1.0)
            else:
                # Default blue to red
                lut.SetHueRange(0.667, 0.0)
                lut.SetSaturationRange(1.0, 1.0)
                lut.SetValueRange(0.3, 1.0)
                
            lut.Build()
            return lut
            
        except Exception as e:
            print(f"Error creating lookup table: {e}")
            # Return simple blue-to-red as fallback
            lut.SetHueRange(0.667, 0.0)
            lut.Build()
            return lut

    def update_opacity(self, value):
        """Update field visualization opacity - FIXED"""
        opacity = value / 100.0
        self.opacity_label.setText(f"{value}%")
        
        try:
            if hasattr(self, 'field_actor') and self.field_actor:
                self.field_actor.GetProperty().SetOpacity(opacity)
                
            if hasattr(self, 'volume_actor') and self.volume_actor:
                # For volume rendering, update the opacity function
                volume_property = self.volume_actor.GetProperty()
                opacity_func = volume_property.GetScalarOpacity()
                
                # Scale the existing opacity function
                if opacity_func:
                    # Get current range
                    range_val = opacity_func.GetRange()
                    if range_val[1] > range_val[0]:
                        # Create new scaled opacity function
                        new_opacity = vtk.vtkPiecewiseFunction()
                        for i in range(int(opacity_func.GetSize())):
                            x = range_val[0] + i * (range_val[1] - range_val[0]) / (opacity_func.GetSize() - 1)
                            y = opacity_func.GetValue(x) * opacity
                            new_opacity.AddPoint(x, y)
                        volume_property.SetScalarOpacity(new_opacity)
                
            if hasattr(self, 'slice_actors'):
                for actor in self.slice_actors:
                    if actor:
                        actor.GetProperty().SetOpacity(opacity)
                        
            self.vtk_widget.GetRenderWindow().Render()
            
        except Exception as e:
            print(f"Error updating opacity: {e}")

    def toggle_threshold(self, enabled):
        """Toggle threshold filtering - FIXED"""
        try:
            if enabled and self.vtk_data:
                self.update_threshold()
            else:
                # Restore original data visualization
                current_mode = self.viz_mode_combo.currentText()
                self.change_visualization_mode(current_mode)
        except Exception as e:
            print(f"Error toggling threshold: {e}")

    def update_threshold(self):
        """Update threshold values - FIXED"""
        if not self.threshold_enabled.isChecked() or not self.vtk_data:
            return
            
        try:
            scalar_range = self.vtk_data.GetScalarRange()
            if scalar_range[1] <= scalar_range[0]:
                print("Invalid scalar range for thresholding")
                return
                
            min_val = scalar_range[0] + (scalar_range[1] - scalar_range[0]) * (self.threshold_min_slider.value() / 100.0)
            max_val = scalar_range[0] + (scalar_range[1] - scalar_range[0]) * (self.threshold_max_slider.value() / 100.0)
            
            if min_val >= max_val:
                print("Invalid threshold range")
                return
                
            print(f"Applying threshold: {min_val:.2e} to {max_val:.2e}")
            
            # Apply threshold filter
            threshold = vtk.vtkThreshold()
            threshold.SetInputData(self.vtk_data)
            threshold.SetLowerThreshold(min_val)
            threshold.SetUpperThreshold(max_val)
            threshold.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)
            threshold.Update()
            
            thresholded_output = threshold.GetOutput()
            print(f"Threshold result: {thresholded_output.GetNumberOfPoints()} points")
            
            if thresholded_output.GetNumberOfPoints() > 0:
                # Apply current visualization mode to thresholded data
                self.apply_visualization_to_data(thresholded_output)
            else:
                print("No points passed threshold")
                
        except Exception as e:
            print(f"Error updating threshold: {e}")

    def apply_visualization_to_data(self, data):
        """Apply current visualization mode to given data - FIXED"""
        if not data:
            return
            
        try:
            # Clear existing visualization
            self.clear_field_visualization()
            
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
                self.setup_field_visualization()
                
            # Restore original data reference
            self.vtk_data = original_data
            
            self.vtk_widget.GetRenderWindow().Render()
            
        except Exception as e:
            print(f"Error applying visualization to data: {e}")
            # Restore original data and try basic visualization
            self.vtk_data = original_data
            self.setup_field_visualization()

    def toggle_isosurfaces(self, enabled):
        """Toggle isosurface display - FIXED"""
        try:
            if enabled and self.vtk_data:
                if self.viz_mode_combo.currentText() != "Isosurfaces":
                    self.viz_mode_combo.setCurrentText("Isosurfaces")
                else:
                    self.setup_isosurface_rendering()
            elif not enabled:
                # Switch to a different mode or clear isosurfaces
                if self.viz_mode_combo.currentText() == "Isosurfaces":
                    self.viz_mode_combo.setCurrentText("Point Cloud")
        except Exception as e:
            print(f"Error toggling isosurfaces: {e}")

    def update_isosurfaces(self):
        """Update number of isosurfaces - FIXED"""
        try:
            if (self.viz_mode_combo.currentText() == "Isosurfaces" and 
                self.isosurface_enabled.isChecked() and 
                self.vtk_data):
                self.setup_isosurface_rendering()
        except Exception as e:
            print(f"Error updating isosurfaces: {e}")

    def change_colormap(self, colormap_name):
        """Change the color mapping - FIXED"""
        try:
            new_lut = self.create_lookup_table(colormap_name)
            
            if hasattr(self, 'field_actor') and self.field_actor:
                mapper = self.field_actor.GetMapper()
                if mapper:
                    mapper.SetLookupTable(new_lut)
                    
            if hasattr(self, 'slice_actors'):
                for actor in self.slice_actors:
                    if actor:
                        mapper = actor.GetMapper()
                        if mapper:
                            mapper.SetLookupTable(new_lut)
                            
            # Update scalar bar
            if hasattr(self, 'scalar_bar') and self.scalar_bar:
                self.scalar_bar.SetLookupTable(new_lut)
                    
            self.vtk_widget.GetRenderWindow().Render()
            print(f"Changed colormap to: {colormap_name}")
            
        except Exception as e:
            print(f"Error changing colormap: {e}")
        
    # def create_earth_representation(self):
    #     """Create a simple Earth sphere for reference"""
    #     earth_sphere = vtk.vtkSphereSource()
    #     earth_sphere.SetRadius(6371.0)  # Earth radius in km
    #     earth_sphere.SetThetaResolution(50)
    #     earth_sphere.SetPhiResolution(50)
        
    #     earth_mapper = vtk.vtkPolyDataMapper()
    #     earth_mapper.SetInputConnection(earth_sphere.GetOutputPort())
        
    #     self.earth_actor = vtk.vtkActor()
    #     self.earth_actor.SetMapper(earth_mapper)
    #     self.earth_actor.GetProperty().SetColor(0.3, 0.3, 0.8)  # Blue-ish
    #     self.earth_actor.GetProperty().SetOpacity(0.3)
        
    #     self.renderer.AddActor(self.earth_actor)
        
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
        """Load VTK data file with comprehensive format support - CLEAN VERSION"""
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
            
            # Use original data directly without conversion for structured grids
            print(f"Using original {type(output).__name__} directly")
            self.vtk_data = output
            
            # Check and setup scalar data
            scalar_array = output.GetPointData().GetScalars()
            if scalar_array:
                print(f"Scalar data: {scalar_array.GetName()}, range: {scalar_array.GetRange()}")
            else:
                print("No primary scalars found, checking available arrays...")
                point_data = output.GetPointData()
                for i in range(point_data.GetNumberOfArrays()):
                    array = point_data.GetArray(i)
                    print(f"  Array {i}: {array.GetName()}")
                
                # Set first array as scalars if available
                if point_data.GetNumberOfArrays() > 0:
                    first_array = point_data.GetArray(0)
                    point_data.SetScalars(first_array)
                    print(f"Set '{first_array.GetName()}' as primary scalars")
                    scalar_array = first_array
            
            # Verify scalar data exists
            if not scalar_array:
                raise ValueError("No scalar data available for visualization")
            
            # Clear any existing visualization
            self.clear_field_visualization()
            
            # Setup field visualization
            print("Setting up field visualization...")
            self.setup_field_visualization()
            
            # Update analyzer
            self.flux_analyzer.set_vtk_data(self.vtk_data)
            
            self.progress_bar.setValue(100)
            
            # Update status
            scalar_name = scalar_array.GetName() if scalar_array else "None"
            scalar_range = scalar_array.GetRange() if scalar_array else (0, 0)
            
            self.status_label.setText(
                f"✓ Loaded VTK data successfully!\n"
                f"Points: {self.vtk_data.GetNumberOfPoints():,} | "
                f"Cells: {self.vtk_data.GetNumberOfCells():,}\n"
                f"Scalar: {scalar_name} | "
                f"Range: {scalar_range[0]:.2e} to {scalar_range[1]:.2e}"
            )
            
            # Force render
            print("Forcing render...")
            self.vtk_widget.GetRenderWindow().Render()
            
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
        """Setup field visualization with automatic zoom to fit data"""
        if not self.vtk_data:
            print("ERROR: No VTK data for field visualization")
            return
            
        print("Setting up field visualization...")
        
        # Determine if we should use point cloud as default
        current_mode = self.viz_mode_combo.currentText() if hasattr(self, 'viz_mode_combo') else "Point Cloud"
        
        if current_mode == "Point Cloud":
            self.setup_point_cloud_rendering()
        else:
            # Use the existing field setup for other modes
            try:
                scalar_array = self.vtk_data.GetPointData().GetScalars()
                if not scalar_array:
                    print("ERROR: No scalar data for visualization")
                    return
                    
                scalar_range = scalar_array.GetRange()
                print(f"Scalar range: {scalar_range[0]:.2e} to {scalar_range[1]:.2e}")
                
                num_points = self.vtk_data.GetNumberOfPoints()
                print(f"Total points: {num_points}")
                
                # Create a threshold filter to only show points with significant flux
                threshold = vtk.vtkThreshold()
                threshold.SetInputData(self.vtk_data)
                threshold.SetLowerThreshold(scalar_range[1] * 0.01)
                threshold.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_UPPER)
                threshold.Update()
                
                thresholded_data = threshold.GetOutput()
                print(f"Points with significant flux: {thresholded_data.GetNumberOfPoints()}")
                
                if thresholded_data.GetNumberOfPoints() == 0:
                    print("No points above threshold, using all points")
                    thresholded_data = self.vtk_data
                
                # Create mapper
                field_mapper = vtk.vtkDataSetMapper()
                field_mapper.SetInputData(thresholded_data)
                field_mapper.SetScalarRange(scalar_range)
                field_mapper.ScalarVisibilityOn()
                
                # Setup lookup table
                lut = self.create_lookup_table(self.colormap_combo.currentText())
                field_mapper.SetLookupTable(lut)
                
                # Create field actor
                self.field_actor = vtk.vtkActor()
                self.field_actor.SetMapper(field_mapper)
                self.field_actor.GetProperty().SetOpacity(0.7)
                
                self.renderer.AddActor(self.field_actor)
                
                # Setup scalar bar
                self.setup_scalar_bar(lut, scalar_array.GetName())
                
                print("Field visualization setup complete")
                
            except Exception as e:
                print(f"ERROR in setup_field_visualization: {e}")
                import traceback
                traceback.print_exc()
        
        # AUTOMATIC ZOOM TO FIT DATA
        self.zoom_to_fit_data()

    def zoom_to_fit_data(self):
        """Automatically zoom to encompass the entire VTK field"""
        if not self.vtk_data:
            return
            
        print("Auto-zooming to fit entire VTK field...")
        
        try:
            # Get data bounds
            bounds = self.vtk_data.GetBounds()
            print(f"Data bounds: X({bounds[0]:.0f}, {bounds[1]:.0f}) Y({bounds[2]:.0f}, {bounds[3]:.0f}) Z({bounds[4]:.0f}, {bounds[5]:.0f})")
            
            # Calculate the center and size of the data
            center = [
                (bounds[1] + bounds[0]) / 2,
                (bounds[3] + bounds[2]) / 2,
                (bounds[5] + bounds[4]) / 2
            ]
            
            # Calculate the maximum extent
            x_range = bounds[1] - bounds[0]
            y_range = bounds[3] - bounds[2]
            z_range = bounds[5] - bounds[4]
            max_range = max(x_range, y_range, z_range)
            
            print(f"Data center: ({center[0]:.0f}, {center[1]:.0f}, {center[2]:.0f})")
            print(f"Max range: {max_range:.0f} km")
            
            # Position camera to see entire data with some margin
            camera = self.renderer.GetActiveCamera()
            
            # Set focal point to data center
            camera.SetFocalPoint(center[0], center[1], center[2])
            
            # Position camera at a distance that shows all data with 20% margin
            distance = max_range * 1.5  # 1.5x for good margin
            
            # Position camera at 45-degree angle for good 3D view
            camera.SetPosition(
                center[0] + distance * 0.7,  # X offset
                center[1] + distance * 0.7,  # Y offset  
                center[2] + distance * 0.5   # Z offset (slightly above)
            )
            
            # Set up vector (Z-axis up)
            camera.SetViewUp(0, 0, 1)
            
            # Reset camera to ensure proper clipping planes
            self.renderer.ResetCamera()
            
            # Fine-tune the zoom to show everything with margin
            camera.Zoom(0.8)  # Zoom out a bit more for safety margin
            
            # Force render
            self.vtk_widget.GetRenderWindow().Render()
            
            print(f"Camera positioned at: {camera.GetPosition()}")
            print("Auto-zoom complete - entire VTK field should be visible")
            
        except Exception as e:
            print(f"Error in auto-zoom: {e}")
            # Fallback to default camera position
            camera = self.renderer.GetActiveCamera()
            camera.SetPosition(20000, 20000, 10000)
            camera.SetFocalPoint(0, 0, 0)
            camera.SetViewUp(0, 0, 1)
            self.renderer.ResetCamera()
        
    def debug_field_actor(self):
        """Specific debug for field actor"""
        print("\n=== FIELD ACTOR DEBUG ===")
        
        if not hasattr(self, 'field_actor') or not self.field_actor:
            print("No field_actor exists!")
            return
            
        print("Field actor exists")
        
        # Check mapper
        mapper = self.field_actor.GetMapper()
        if mapper:
            print("Mapper exists")
            input_data = mapper.GetInput()
            if input_data:
                print(f"Mapper input: {input_data.GetNumberOfPoints()} points")
                scalar_range = mapper.GetScalarRange()
                print(f"Mapper scalar range: {scalar_range}")
            else:
                print("No input data on mapper!")
        else:
            print("No mapper on field actor!")
            
        # Check properties
        prop = self.field_actor.GetProperty()
        print(f"Opacity: {prop.GetOpacity()}")
        print(f"Point size: {prop.GetPointSize()}")
        print(f"Color: {prop.GetColor()}")
        print(f"Representation: {prop.GetRepresentation()}")
        
        # Check visibility
        print(f"Visibility: {self.field_actor.GetVisibility()}")
        
        print("=========================\n")
            
    def setup_scalar_bar(self, lut, scalar_name):
        """ENHANCED: Setup scalar bar with proper units"""
        try:
            # Remove existing scalar bar
            if hasattr(self, 'scalar_bar') and self.scalar_bar:
                self.renderer.RemoveViewProp(self.scalar_bar)
                
            self.scalar_bar = vtk.vtkScalarBarActor()
            self.scalar_bar.SetLookupTable(lut)
            
            # Enhanced title with units
            if scalar_name:
                if "flux" in scalar_name.lower():
                    title = "Electron Flux\n(particles/cm²/s)"
                else:
                    title = f"{scalar_name}\n(data units)"
            else:
                title = "Field Value\n(data units)"
                
            self.scalar_bar.SetTitle(title)
            self.scalar_bar.SetPosition(0.85, 0.1)
            self.scalar_bar.SetWidth(0.12)
            self.scalar_bar.SetHeight(0.8)
            
            # Improve scalar bar appearance
            self.scalar_bar.SetNumberOfLabels(6)
            self.scalar_bar.GetLabelTextProperty().SetColor(1, 1, 1)
            self.scalar_bar.GetTitleTextProperty().SetColor(1, 1, 1)
            self.scalar_bar.GetTitleTextProperty().SetFontSize(11)
            self.scalar_bar.GetLabelTextProperty().SetFontSize(9)
            
            # Format numbers in scientific notation for better readability
            self.scalar_bar.GetLabelTextProperty().SetJustificationToLeft()
            
            self.renderer.AddViewProp(self.scalar_bar)
            print("Enhanced scalar bar with units added")
            
        except Exception as e:
            print(f"Error setting up enhanced scalar bar: {e}")

    def update_status_with_units(self):
        """ENHANCED: Update status display with proper units"""
        if not self.orbital_path or self.current_time_index >= len(self.orbital_path):
            return
            
        current_point = self.orbital_path[self.current_time_index]
        
        # Calculate additional useful metrics
        distance_from_earth = np.sqrt(current_point.x**2 + current_point.y**2 + current_point.z**2)
        altitude_km = distance_from_earth - 6371  # Earth radius
        
        # Calculate flux with units
        flux = self.flux_analyzer.analyze_flux_at_point(current_point)
        
        # Enhanced status with better units and formatting
        self.status_label.setText(
            f"Time: {current_point.time:.2f} h | Altitude: {altitude_km:.1f} km\n"
            f"Position: ({current_point.x:.0f}, {current_point.y:.0f}, {current_point.z:.0f}) km\n"
            f"Flux: {flux:.2e} particles/s through {self.cross_section_spinbox.value():.1f} m² cross-section"
        )
        
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
