# Import necessary libraries
import os
import math
import ufl
import matplotlib.pyplot as plt
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, io, log
from dolfinx.common import Timer
from dolfinx.io import VTXWriter
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc
import numpy as np
from excitations import SupplyCurrentDensity, PMMagnetization

# Output directory
output_dir = "./output"
os.makedirs(output_dir, exist_ok=True)

# Helper methods
def AssembleSystem(a, L, bcs, name, mesh):
    """Assemble the system matrix and vector"""
    if MPI.COMM_WORLD.rank == 0:
        log.log(log.LogLevel.WARNING, name + ": Assembling LHS Matrix")
    A = assemble_matrix(a, bcs)
    A.assemble()

    if MPI.COMM_WORLD.rank == 0:
        log.log(log.LogLevel.WARNING, name + ": Assembling RHS Vector")
    b = assemble_vector(L)
    apply_lifting(b, [a], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, bcs)
    return A, b

def SaveSolution(fnc, t, writer):
    """Save the solution for post-processing"""
    fnc.x.scatter_forward()
    writer.write(t)

# Derived quantities for post-processing
class DerivedQuantitiesPMSM:
    def __init__(self, A, mesh, domains, sigma, omega, x):
        self.A = A
        self.mesh = mesh
        self.domains = domains
        self.sigma = sigma
        self.omega = omega
        self.x = x  # Spatial coordinates
        self.mu_0 = 1.25663706143592e-06

    def compute_torque(self, dt):
        """
        Compute torque using Maxwell stress tensor approach.
        """
        r = ufl.sqrt(self.x[0]**2 + self.x[1]**2)
        torque_form = (
            self.sigma
            * self.omega
            * ufl.dot(ufl.as_vector((-self.x[1], self.x[0])), ufl.grad(self.A))
            * r
            * dx(self.domains["Rotor"])
        )
        torque = fem.assemble_scalar(fem.form(torque_form))
        return self.mesh.comm.allreduce(torque, op=MPI.SUM)

# Parameters
freq = 50.0
motorrpm = 600.0
omega_v = motorrpm / 9.5492965964254
omega_f = 2.0 * math.pi * freq
mu_0 = 1.25663706143592e-06
sigma_stl = 2.00e+06

t = 0.0
dt = 0.001
t_final = 0.005

# Mesh and function spaces
meshname = "meshes/test4.xdmf"
with io.XDMFFile(MPI.COMM_WORLD, meshname, "r") as file:
    mesh = file.read_mesh()
    mt_domains = file.read_meshtags(mesh, "Cell_markers")

dx = ufl.Measure("dx", subdomain_data=mt_domains, domain=mesh)
V_A = fem.FunctionSpace(mesh, ("CG", 1))
V_V = fem.FunctionSpace(mesh, ("CG", 1))
V_VF = fem.FunctionSpace(mesh, ufl.MixedElement([V_A.ufl_element(), V_V.ufl_element()]))

# Initialize functions
u_a, u_v = ufl.TrialFunctions(V_VF)
v_a, v_v = ufl.TestFunctions(V_VF)
A0 = fem.Function(V_A)
A_sol = fem.Function(V_A, name="Az")
V_sol = fem.Function(V_V, name="V")
Az_vtx = VTXWriter(mesh.comm, os.path.join(output_dir, "A.bp"), [A_sol], engine="BP4")
Vz_vtx = VTXWriter(mesh.comm, os.path.join(output_dir, "V.bp"), [V_sol], engine="BP4")

# Excitations
jsexp = SupplyCurrentDensity()
jsexp.amp = 1.0  # Set current amplitude
jsource = fem.Function(V_A)
jsource.interpolate(jsexp.eval)

msexp = PMMagnetization()
msexp.mag = 1.0
msource = fem.Function(V_A)
msource.interpolate(msexp.eval)

# Variational forms
f_a = (
    (1 / mu_0) * ufl.inner(ufl.grad(u_a), ufl.grad(v_a)) * dx
    + sigma_stl * ufl.dot(ufl.grad(u_v), ufl.grad(v_a)) * dx
    - ufl.dot(jsource, ufl.grad(v_a)) * dx
    - ufl.inner(msexp.mag * ufl.as_vector((u_a.dx(1), -u_a.dx(0))), v_a) * dx
)
f_v = (
    sigma_stl * ufl.inner(ufl.grad(u_v), ufl.grad(v_v)) * dx
    - (1 / dt) * sigma_stl * (u_a - A0) * ufl.div(ufl.grad(v_v)) * dx
)
form_av = f_a + f_v
a_av, L_av = ufl.system(form_av)

# Assemble system
bcs_av = []
A_form = fem.form(a_av)
L_form = fem.form(L_av)
A_av = assemble_matrix(A_form, bcs_av)
solver_av = PETSc.KSP().create(mesh.comm)
solver_av.setOperators(A_av)
solver_av.setType("gmres")
solver_av.getPC().setType("bjacobi")
solver_av.setFromOptions()

# Derived quantities
x = ufl.SpatialCoordinate(mesh)
domains = {"Rotor": (5,)}  # Replace with your domain IDs
derived = DerivedQuantitiesPMSM(A_sol, mesh, domains, sigma_stl, omega_v, x)

# Time-stepping
torques = []
times = []
while t <= t_final:
    # Assemble system
    A_av, b_av = AssembleSystem(A_form, L_form, bcs_av, "AV", mesh)

    # Solve
    solver_av.solve(b_av, A_sol.vector)
    A_sol.x.scatter_forward()
    V_sol.x.scatter_forward()

    # Save solutions
    SaveSolution(A_sol, t, Az_vtx)
    SaveSolution(V_sol, t, Vz_vtx)

    # Compute torque
    torque = derived.compute_torque(dt)
    torques.append(torque)
    times.append(t)

    # Update
    A0.x.array[:] = A_sol.x.array[:]
    t += dt

# Close writers
Az_vtx.close()
Vz_vtx.close()

# Plot torque vs time
if MPI.COMM_WORLD.rank == 0:
    plt.figure()
    plt.plot(times, torques, label="Torque")
    plt.xlabel("Time (s)")
    plt.ylabel("Torque (N·m)")
    plt.title("Torque vs Time")
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(output_dir, "torque_vs_time.png"))
    plt.show()
