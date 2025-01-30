# from meshinfo import PmsmMeshValues
# from excitations import SupplyCurrentDensity, PMMagnetization

import os, math, ufl
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
from io import TextIOWrapper
from dolfinx import cpp, fem, io, default_scalar_type

from ufl import (FiniteElement, VectorElement, MixedElement, 
    TestFunction, TrialFunction, TestFunctions, TrialFunctions,
    SpatialCoordinate, grad, cross, curl, inner, system, 
    as_vector, FacetNormal, dot, sqrt, avg)

from dolfinx import log
# from dolfinx import (Function, FunctionSpace, log,
    # TimingType, list_timings)       # DirichletBC
from dolfinx.common import Timer
from dolfinx.cpp.mesh import GhostMode
from dolfinx.io import XDMFFile
from dolfinx.fem import (assemble_matrix, assemble_vector,
    assemble_scalar, apply_lifting, set_bc, locate_dofs_topological)

t_run = Timer("00 Overall Run Time")

# Log output options
# loglevel = log.LogLevel.INFO
loglevel = log.LogLevel.WARNING
# loglevel = log.LogLevel.ERROR
log.set_log_level(loglevel)
log.set_output_file("output.txt")

### Parameters ###########################################################

t_01 = Timer("01 Initialise Variables")
if MPI.COMM_WORLD.rank == 0:
    log.log(loglevel, "Defining Model Parameters")

### Add Mesh Name
meshname = "test4.xdmf"

# Frequency (Hz)
freq = 50.0

# Rotational velocity (RPM)
motorrpm = 600.0

# Rotational velocity (rad/s)
omega_v = motorrpm/9.5492965964254

# Angular frequency (Hz)
omega_f = 2.*math.pi*freq

# Relative permeability (-)
mur_air = 1.0
mur_cop = 0.999991
mur_mag = 1.04457
mur_stl = 100.0

# Permeability of free space (H/m)
mu_0 = 1.25663706143592e-06

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

# Supply current (A)
sp_current  = 35.00
# Current density magnitude (A/m^2)
jsource_amp = sp_current/2.47558E-05

# Remanent magnetic flux density (T)
msource_mag_T = 1.09999682447133
# Permanent Magnetization (A/m)
msource_mag   = (msource_mag_T*1e7)/(4*math.pi)

# Uncomment to neglect respective source term     
# jsource_amp = 1.00e-32
# msource_mag = 1.00e-32

# Initial time / Time step / Final time (s)
t         = 0.000
dt        = 0.001
t_final   = 0.005
# dt      = 0.2 if freq == 0 else 0.02/freq
# t_final = 0.8 if freq == 0 else (0.02/freq)*50

# Quadrature degree (-)
quaddeg = 3

# Vector / Scalar order (-)
order_v = 1
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

# with XDMFFile(MPI.COMM_WORLD,
#               meshpath,
#               "r",
#               encoding=XDMFFile.Encoding.HDF5) as file:
#     mesh = file.read_mesh(ghost_mode=GhostMode.none,
#                           name="mesh",
#                           xpath=r"/Xdmf/Domain")
#     # mesh.topology.create_connectivity() # create_connectivity_all()
#     tdim = 2
#     mesh.topology.create_connectivity(tdim - 1, 0)
#     # ct = xdmf.read_meshtags(mesh, name="Cell_markers")
#     mt_facets  = file.read_meshtags(mesh, "facets")
#     mt_domains = file.read_meshtags(mesh, "domains")

# Read mesh and cell markers
with io.XDMFFile(MPI.COMM_WORLD, meshpath, "r") as xdmf:
    mesh = xdmf.read_mesh()     # <dolfinx.mesh.Mesh object at 0x7f556028faf0>
    ct = xdmf.read_meshtags(mesh, name="Cell_markers")
    tdim = mesh.topology.dim        # 2
    mesh.topology.create_connectivity(tdim - 1, 0)
    ft = xdmf.read_meshtags(mesh, name="Facet_markers")


# dx = ufl.Measure("dx", subdomain_data=mt_domains, domain=mesh) \
#                 (metadata={"quadrature_degree": quaddeg})
ds = ufl.Measure("ds", subdomain_data=ft,  domain=mesh) \
                (metadata={"quadrature_degree": quaddeg})
# dS = ufl.Measure("dS", subdomain_data=mt_facets,  domain=mesh) \
#                 (metadata={"quadrature_degree": quaddeg})

# mshval = PmsmMeshValues(meshname)
# OutputMesh(mesh, mt_domains, mt_facets, "mesh22")
t_02.stop()