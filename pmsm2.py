import argparse
import sys
import os
from io import TextIOWrapper
from pathlib import Path
from typing import Callable, Optional, TextIO, Union, Dict

import dolfinx.fem.petsc as _petsc
import dolfinx.mesh
import matplotlib.pyplot as plt
import numpy as np
import math
import tqdm
import ufl
from dolfinx import cpp, log, fem, io, default_scalar_type
from dolfinx.io import VTXWriter
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx.common import Timer

from generate_pmsm_2D import (domain_parameters, model_parameters,
                                    surface_map)
from utils2 import DerivedQuantities2D, MagneticField2D, update_current_density
    
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
    if mesh.mpi_comm().rank == 0: 
        log.log(log.LogLevel.WARNING, name + ": Assembling LHS Matrix")
    A = _petsc.assemble_matrix(a, bcs)
    A.assemble()

    if mesh.mpi_comm().rank == 0:
        log.log(log.LogLevel.WARNING, name + ": Assembling RHS Vector")
    b = _petsc.assemble_vector(L)
    _petsc.apply_lifting(b, [a], [bcs])
    b.scatter_reverse(PETSc.ScatterMode.ADD)
    fem.set_bc(b, bcs)
    return A, b

def SaveSolution(fnc, t, file, name, units):
    fnc.x.scatter_forward()
    file.write_function(fnc, round(t,3))
    fncmax, fncmin = fnc.x.array.max(), fnc.x.array.min()
    if fnc.function_space.mesh.mpi_comm().rank == 0:
        log.log(log.LogLevel.WARNING, name + " Max: {:.3e} {}".format(fncmax, units))
        log.log(log.LogLevel.WARNING, name + " Min: {:.3e} {}".format(fncmin, units))



t_run = Timer("00 Overall Run Time")

# Log output options
# loglevel = log.LogLevel.INFO
loglevel = log.LogLevel.WARNING
# loglevel = log.LogLevel.ERROR
log.set_log_level(loglevel)
log.set_output_file("output.txt")

t_01 = Timer("01 Initialise Variables")
if MPI.COMM_WORLD.rank == 0:
    log.log(loglevel, "Defining Model Parameters")

### Parameters ###########################################################

meshname = "test4.xdmf"                 # Mesh File Name
freq = 50.0                             # Frequency (Hz)
motorrpm = 600.0                        # Rotational velocity (RPM)
omega_v = motorrpm/9.5492965964254      # Rotational velocity (rad/s)
omega_f = 2.*math.pi*freq               # Angular frequency (Hz)

mu_0 = 1.25663706143592e-06             # Permeability of free space (H/m)

# Relative permeability (-)
mur_air = 1.0
mur_cop = 0.999991
mur_mag = 1.04457
mur_stl = 100.0

# Permeability (H/m)
mu_air = mur_air * mu_0
mu_cop = mur_cop * mu_0
mu_mag = mur_mag * mu_0
mu_stl = mur_stl * mu_0

# Electrical conductiviy (S/m)
sigma_air = 1.00e-32
sigma_cop = 1.00e-32
sigma_mag = 1.00e-32
sigma_stl = 2.00e+06


sp_current  = 35.00                                 # Supply current (A)
jsource_amp = sp_current/2.47558E-05                # Current density magnitude (A/m^2
msource_mag_T = 1.09999682447133                    # Remanent magnetic flux density (T)
msource_mag   = (msource_mag_T*1e7)/(4*math.pi)     # Permanent Magnetization (A/m)


t         = 0.000           # Initial time (s)
dt        = 0.001           # Time step (s)
t_final   = 0.005           # Final time (s)

# dt      = 0.2 if freq == 0 else 0.02/freq
# t_final = 0.8 if freq == 0 else (0.02/freq)*50


quaddeg = 3             # Quadrature degree (-)
order_v = 1             # Vector / Scalar order (-)
order_s = 1

# Helper params
resdirname   = "Results_" + str(int(freq)) + "Hz"
skiptimeloop = False

t_01.stop()

### Mesh #################################################################

t_02 = Timer("02 Read Mesh")
filedir = os.getcwd()
meshpath = os.path.join(filedir, "meshes", meshname)
if MPI.COMM_WORLD.rank == 0: 
    log.log(loglevel, "Reading Mesh: " + meshpath)

# Read mesh and cell markers
with io.XDMFFile(MPI.COMM_WORLD, meshpath, "r") as xdmf:
    mesh = xdmf.read_mesh()     # <dolfinx.mesh.Mesh object at 0x7f556028faf0>
    ct = xdmf.read_meshtags(mesh, name="Cell_markers")
    tdim = mesh.topology.dim        # 2
    mesh.topology.create_connectivity(tdim - 1, 0)
    ft = xdmf.read_meshtags(mesh, name="Facet_markers")

# with XDMFFile(MPI.COMM_WORLD,
#               meshpath,
#               "r",
#               encoding=XDMFFile.Encoding.HDF5) as file:
#     mesh = file.read_mesh(ghost_mode=GhostMode.none,
#                           name="mesh",
#                           xpath=r"/Xdmf/Domain")
#     mesh.topology.create_connectivity_all()
#     mt_facets  = file.read_meshtags(mesh, "facets")
#     mt_domains = file.read_meshtags(mesh, "domains")

# dx = ufl.Measure("dx", subdomain_data=ct, domain=mesh) \
#                 (metadata={"quadrature_degree": quaddeg})       # ct = mt_domains
# ds = ufl.Measure("ds", subdomain_data=ft,  domain=mesh) \
#                 (metadata={"quadrature_degree": quaddeg})
# dS = ufl.Measure("dS", subdomain_data=ft,  domain=mesh) \
#                 (metadata={"quadrature_degree": quaddeg})

domains: Dict[str, tuple[int, ...]] = {"Cu": (7, 8, 9, 10, 11, 12), "Stator": (6, ), "Rotor": (5, ),
                                                 "Al": (4,), "AirGap": (2, 3), "Air": (1,), "PM": (13, 14, 15, 16, 17, 18, 19, 20, 21, 22), "Al2": (23,)}

# mshval = PmsmMeshValues(meshname)
# OutputMesh(mesh, mt_domains, mt_facets, "mesh22")
t_02.stop()

### Function Spaces ######################################################

if MPI.COMM_WORLD.rank == 0: log.log(loglevel, "Creating Function Spaces")
# Define problem function space
cell = mesh.ufl_cell()  # triangle
degree = 1
FE = ufl.FiniteElement("Lagrange", cell, degree)
ME = ufl.MixedElement([FE, FE])
VQ = fem.FunctionSpace(mesh, ME)  # VQ: This is the mixed function space containing both the magnetic vector potential component A_z and the scalar potential V.

# Define test, trial and functions for previous timestep
Az, V = ufl.TrialFunctions(VQ)
vz, q = ufl.TestFunctions(VQ)
AnVn = fem.Function(VQ)
An, _ = ufl.split(AnVn)  # Solution at previous time step
DG0 = fem.FunctionSpace(mesh, ("DG", 0))
J0z = fem.Function(DG0)  # Current density
## Add magnetic component

# Create integration sets
# Omega = domains["Cu"] + domains["Stator"] + domains["Air"] + domains["AirGap"]
# Omega_c = domains["Rotor"] + domains["Al"]

dx_air = domains["Air"] + domains["AirGap"]
dx_rotor = domains["Rotor"] + domains["Al"] + domains["Al2"]
dx_stator = domains["Stator"]
dx_coil = domains["Cu"]
dx_magnet = domains["PM"]
dx_all = dx_air + dx_rotor + dx_stator + dx_coil + dx_magnet
# Create integration measures
dx = ufl.Measure("dx", domain=mesh, subdomain_data=ct)
ds = ufl.Measure("ds", domain=mesh, subdomain_data=ft)  # subdomain_id=surface_map["Exterior"])


# Define temporal and spatial parameters
n = ufl.FacetNormal(mesh)
dt = fem.Constant(mesh, dt)
x = ufl.SpatialCoordinate(mesh)

omega = fem.Constant(mesh, default_scalar_type(omega_v))

# Define variational form
a = dt / mu_0 * ufl.inner(ufl.grad(Az), ufl.grad(vz)) * dx(dx_all)          # check mu values
a += mu_0 * sigma * Az * vz * dx(dx_rotor + dx_stator)


##########


