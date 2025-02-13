'''
Python script to convert .bp to .vtu to analyze in Paraview
'''

import adios2
import numpy as np
import pyvista as pv

# Open the .bp file
f_loc = "PMSM_Jan_09_11_34_44/B"

with adios2.open(f"{f_loc}.bp", "r") as fh:
    # Extract data from the file
    for step in fh:
        # Mesh data
        geometry = step.read("geometry").reshape(-1, 3)  # Reshape to (N, 3) for 3D coordinates
        connectivity = step.read("connectivity")  # Connectivity array
        types = step.read("types")  # VTK cell types
        az_data = step.read("B")  # Field data Az, B

if np.unique(types).size > 1:
    raise ValueError("Multiple cell types detected. This script handles one type at a time.")

# Build the mesh using PyVista
# PyVista expects the connectivity to start with the number of points per cell
# Format: [n_points, p1, p2, ..., n_points, p1, p2, ...]
connectivity_fixed = connectivity[:, 1:]


# Add cell size as the first value of each element for PyVista
n_cells = connectivity_fixed.shape[0]
cells = np.insert(connectivity_fixed, 0, 3, axis=1).flatten()  # Prepend "3" for each cell

# Create an array of cell types (one per cell)
cell_types = np.full(n_cells, types, dtype=np.uint8)

# Create an unstructured grid
grid = pv.UnstructuredGrid(cells, cell_types, geometry)

# Add the field data 
grid["B"] = az_data # Az, B

# Save the mesh as a .vtu file
grid.save(f"{f_loc}.vtu")
print(f"Saved to {f_loc}.vtu")
