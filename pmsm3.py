# from meshinfo import PmsmMeshValues
from excitations import SupplyCurrentDensity, PMMagnetization

import os, math, ufl
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
from typing import Dict
from ufl import (FiniteElement, VectorElement, MixedElement, 
                 TestFunction, TrialFunction, TestFunctions, TrialFunctions,
                 SpatialCoordinate, div, grad, curl, inner, system, 
                 as_vector, FacetNormal, dot, sqrt, avg, dx, ds, dS)
from dolfinx import cpp, fem, io, default_scalar_type, log
from dolfinx.common import Timer
from dolfinx.cpp.mesh import GhostMode
from dolfinx.io import XDMFFile
from dolfinx.fem import (assemble_scalar, locate_dofs_topological)
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc
from dolfinx.io import VTXWriter
from dolfinx.common import list_timings, TimingType


output_dir = "./output"
os.makedirs(output_dir, exist_ok=True)



t_run = Timer("00 Overall Run Time")
log.set_log_level(log.LogLevel.WARNING)

### Helper Methods
def SumFunctions(fncspace, name, A, B):
    fnc = fem.Function(fncspace)
    fnc.name = name
    fnc.x.array.set(0.0)
    fnc.x.scatter_forward()
    fnc.x.axpy(1, A.x)
    fnc.x.axpy(1, B.x)
    return fnc

def AssembleSystem(a, L, bcs, name, mesh):
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

def SaveSolution(fnc, t, file, name, units):
    fnc.x.scatter_forward()
    file.write_function(fnc, round(t,3))
    fncmax, fncmin = fnc.x.array.max(), fnc.x.array.min()
    if fnc.function_space.mesh.mpi_comm().rank == 0:
        log.log(log.LogLevel.WARNING, name + " Max: {:.3e} {}".format(fncmax, units))
        log.log(log.LogLevel.WARNING, name + " Min: {:.3e} {}".format(fncmin, units))



### Parameters
t_01 = Timer("01 Initialise Variables")
meshname = "test4.xdmf"  # 2D mesh file
freq = 50.0
motorrpm = 600.0
omega_v = motorrpm/9.5492965964254
omega_f = 2.*math.pi*freq

mu_0 = 1.25663706143592e-06
mur_air = 1.0
mur_stl = 100.0
mur_mag = 1.04457
mur_cop = 0.999991

mu_air = mur_air * mu_0
mu_cop = mur_cop * mu_0
mu_mag = mur_mag * mu_0
mu_stl = mur_stl * mu_0

sigma_air = 1.00e-32
sigma_cop = 1.00e-32
sigma_mag = 1.00e-32
sigma_stl = 2.00e+06

sp_current  = 35.00
jsource_amp = sp_current/2.47558E-05
msource_mag_T = 1.09999682447133
msource_mag   = (msource_mag_T*1e7)/(4*math.pi)

t = 0.000
dt = 0.001
t_final = 0.005
quaddeg = 3
order_v = 1
order_s = 1
t_01.stop()

### Mesh (2D)
t_02 = Timer("02 Read Mesh")
filedir = os.getcwd()
meshpath = os.path.join(filedir, "meshes", meshname)
if MPI.COMM_WORLD.rank == 0: 
    log.log(log.LogLevel.WARNING, "Reading Mesh: " + meshpath)

with XDMFFile(MPI.COMM_WORLD,
              meshpath,
              "r",
              encoding=XDMFFile.Encoding.HDF5) as file:
    mesh = file.read_mesh(ghost_mode=GhostMode.none)
    mesh.topology.create_connectivity(1, 2)
    mt_facets  = file.read_meshtags(mesh, "Facet_markers")
    mt_domains = file.read_meshtags(mesh, "Cell_markers")   # domains -- check--

dx = ufl.Measure("dx", subdomain_data=mt_domains, domain=mesh)(metadata={"quadrature_degree": quaddeg})
ds = ufl.Measure("ds", subdomain_data=mt_facets,  domain=mesh)(metadata={"quadrature_degree": quaddeg})
dS = ufl.Measure("dS", subdomain_data=mt_facets, domain=mesh)(metadata={"quadrature_degree": quaddeg})

# mshval = PmsmMeshValues(meshname)
t_02.stop()

### Function Spaces (2D)
if MPI.COMM_WORLD.rank == 0: 
    log.log(log.LogLevel.WARNING, "Creating Function Spaces")

cell = mesh.ufl_cell()

# In a 2D magnetostatic model, A typically has one out-of-plane component (A_z). Representing A as a scalar:
CG_F  = FiniteElement("CG", cell, order_s)  # For A as scalar
CG_V  = VectorElement("CG", cell, order_v)  # For vector fields if needed
CG_VF = MixedElement([CG_F, CG_F]) # (A,V) both scalars if modeling A_z and a scalar V

V_V    = fem.FunctionSpace(mesh, CG_V)   # For vector fields like J or M if needed
V_VF   = fem.FunctionSpace(mesh, CG_VF)  # For (A, V) both as scalars
V_DGV  = fem.FunctionSpace(mesh, VectorElement("DG", cell, 0))

### UFL Functions
(x,y) = SpatialCoordinate(mesh)
normal = FacetNormal(mesh)
vec0   = as_vector((0,0))
z_unit = as_vector((0,0)) # In 2D, no z dimension, but keep as dummy if needed.

### Trial and Test Functions
t_03 = Timer("03 Define Problem")
u_a, u_v = TrialFunctions(V_VF)
v_a, v_v = TestFunctions(V_VF)

# Create current excitation
DG0 = fem.FunctionSpace(mesh, ("DG", 0))
jsexp = SupplyCurrentDensity() #fem.Function(DG0)   #SupplyCurrentDensity()         --check--
jsexp.amp = jsource_amp
jsexp.omega = omega_f
jsource = fem.Function(V_V)
jsource.interpolate(jsexp.eval)
jsource.x.scatter_forward()

# Create magnetization excitation
msexp = PMMagnetization()   #fem.Function(DG0)   #PMMagnetization()
msexp.mag = msource_mag
msource = fem.Function(V_V)
msource.interpolate(msexp.eval)
msource.x.scatter_forward()

# Initial A
V_A_collapse, injective_map = V_VF.sub(0).collapse()
A0 = fem.Function(V_A_collapse) # V_VF.sub(0).collapse()
A0.x.array[:] = 0.0
A0.x.scatter_forward()

# Domain labeling (as in original code, adapt IDs for 2D mesh)
domains: Dict[str, tuple[int, ...]] = {"Cu": (7, 8, 9, 10, 11, 12), "Stator": (6, ), "Rotor": (5, ),
                                                 "Al": (4,), "AirGap": (2, 3), "Air": (1,), "PM": (13, 14, 15, 16, 17, 18, 19, 20, 21, 22), "Al2": (23,)}
dx_air = domains["Air"] + domains["AirGap"]
dx_rotor = domains["Rotor"] + domains["Al"] + domains["Al2"]
dx_stator = domains["Stator"]
dx_coil = domains["Cu"]
dx_magnet = domains["PM"]
steel_domain = dx(dx_rotor + dx_stator)
rotor_steel_domain = dx(dx_rotor)   
coil_domain  = dx(dx_coil)
magnet_domain= dx(dx_magnet)
air_domain   = dx(dx_air)

# For 2D: curl(A) if A is scalar: curl(A) = (dA/dy, -dA/dx)
# x, y = ufl.SpatialCoordinate(mesh)

def curl_scalar(A):
    return as_vector((A.dx(1), -A.dx(0)))

A = u_a
V = u_v

f_a = inner((1/mu_stl)*grad(A), grad(v_a))*steel_domain \
    + inner((1/mu_cop)*grad(A), grad(v_a))*coil_domain \
    + inner((1/mu_mag)*grad(A), grad(v_a))*magnet_domain \
    + inner((1/mu_air)*grad(A), grad(v_a))*air_domain \
\
    + (1/dt)*inner(sigma_stl*(A - A0), v_a)*steel_domain \
\
    + sigma_stl * dot(grad(V), grad(v_a)) * steel_domain \
    - dot(jsource, grad(v_a)) * coil_domain \
    - inner((mu_0/mu_mag)*msource, curl_scalar(v_a))*magnet_domain \
\
    + sigma_stl * omega_v * dot(as_vector((y, -x)), grad(A)) * v_a * rotor_steel_domain \
    + sigma_mag * omega_v * dot(as_vector((y, -x)), grad(A)) * v_a * magnet_domain


f_v = inner(sigma_stl*grad(V),grad(v_v))*steel_domain \
    + inner(sigma_cop*grad(V),grad(v_v))*coil_domain \
    + inner(sigma_mag*grad(V),grad(v_v))*magnet_domain \
    + inner(sigma_air*grad(V),grad(v_v))*air_domain \
    - (1/dt)*sigma_stl*(A - A0)*div(grad(v_v))*steel_domain \
    + sigma_stl * omega_v * (y * A.dx(0) - x * A.dx(1)) * v_v * rotor_steel_domain \
    + sigma_mag * omega_v * (y * A.dx(0) - x * A.dx(1)) * v_v * magnet_domain


form_av = f_a + f_v
a_av, L_av = system(form_av)

# BCs (2D)
if MPI.COMM_WORLD.rank == 0:
    log.log(log.LogLevel.WARNING, "AV: Applying Boundary Conditions")

bcs_av = []
u0 = fem.Function(V_VF)
u0.x.array[:] = 0.0
u0.x.scatter_forward()

# Suppose we fix A and V=0 on some boundary facet sets:
# Example: If mt_facets has certain IDs for boundaries:
# We'll say facets with ID=1 is where A=0 and V=0
bc_facets = mt_facets.indices[mt_facets.values == 1]
bc_dofs = locate_dofs_topological(V_VF, mesh.topology.dim-1, bc_facets)
bcs_av.append(fem.dirichletbc(u0, bc_dofs))

t_03.stop()

# Setup solver
A_form = fem.form(a_av)
L_form = fem.form(L_av)
A_av = assemble_matrix(A_form, bcs_av)
# A_av.assemble()
# AV = fem.Function(V_VF)

solver_av = PETSc.KSP().create(mesh.comm)
solver_av.setOptionsPrefix("AV_")
solver_av.setOperators(A_av)
solver_av.setTolerances(rtol=1e-08, max_it=100000)
solver_av.setType("gmres")
solver_av.getPC().setType("bjacobi")
solver_av.setFromOptions()

# # Output files

V_collapse, _ = V_VF.sub(0).collapse()

# Create a function in the collapsed subspace
A_sol = fem.Function(V_collapse)
A_sol.name = "A"

# Similarly for V_sol
V_collapse, _ = V_VF.sub(1).collapse()
V_sol = fem.Function(V_collapse)
V_sol.name = "V"

# A_sol = fem.Function(V_VF.sub(0).collapse())
# V_sol = fem.Function(V_VF.sub(1).collapse())
Az_vtx = VTXWriter(mesh.comm, os.path.join(output_dir, "A.bp"), [A_sol], engine="BP4")
Vz_vtx = VTXWriter(mesh.comm, os.path.join(output_dir, "V.bp"), [V_sol], engine="BP4")

t=0.0
while t <= t_final:
    if MPI.COMM_WORLD.rank == 0:
        log.log(log.LogLevel.WARNING, "Solving at t=" + str(t))

    A_av, b_av = AssembleSystem(A_form, L_form, bcs_av, "AV", mesh)

    solver_av.setOperators(A_av)
    x_av = A_av.createVecRight()
    solver_av.solve(b_av, x_av)

    # Convert solution x_av to a DOLFINx Function
    x_av_array = x_av.array  # Access the PETSc Vec data as a NumPy array

    # Split the solution vector into components (A and V)
    # A_sol = fem.Function(V_VF.sub(0).collapse())
    # V_sol = fem.Function(V_VF.sub(1).collapse())

    # Extract data for each component
    # Properly collapse the subspace and extract the function space
    V_collapse_A, _ = V_VF.sub(0).collapse()
    V_collapse_V, _ = V_VF.sub(1).collapse()

    # Access the local DOF information
    dofs_A = V_collapse_A.dofmap.index_map.size_local * V_collapse_A.dofmap.index_map_bs
    dofs_V = V_collapse_V.dofmap.index_map.size_local * V_collapse_V.dofmap.index_map_bs

    A_sol.x.array[:] = x_av_array[:dofs_A]
    V_sol.x.array[:] = x_av_array[dofs_A:]
    
    Az_vtx.write(t)
    Vz_vtx.write(t)

    # Update A0 for the next time step
    A_sol.x.array[:] = A0.x.array[:]  # Copy values
    A0.x.scatter_forward()
    t += dt

# Clean up the VTX writers
Az_vtx.close()
Vz_vtx.close()

if MPI.COMM_WORLD.rank == 0:
    log.log(log.LogLevel.WARNING, "Simulation Finished")

t_run.stop()
list_timings(mesh.comm, [TimingType.wall])


## For visualization see from line 249 
# # Post proc variables of team_30_phi
