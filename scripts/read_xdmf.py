import dolfinx
import meshio
from mpi4py import MPI

# Load the mesh from the XDMF file
with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "mesh1.xdmf", "r") as xdmf_file:
    mesh = xdmf_file.read_mesh()

# Extract cell connectivity for 3D elements
cells = mesh.topology.connectivity(mesh.topology.dim, 0).array
points = mesh.geometry.x

# Convert to Meshio mesh and save as .msh
meshio_mesh = meshio.Mesh(
    points=points,
    cells=[("tetra", cells.reshape(-1, 4))]  # Assuming tetrahedral cells
)
meshio_mesh.write("mesh1.msh")
