"""
Solving a Permanent Magnet Synchronous Motor (PMSM) problem 
using the TEAM-30 code as the foundational framework.
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
from typing import Union, Dict
import math
import time
import dolfinx.fem.petsc as _petsc
import dolfinx.mesh
import numpy as np
import tqdm
import ufl
import basix
from dolfinx import cpp, fem, io, default_scalar_type
from dolfinx.io import VTXWriter
from mpi4py import MPI
from petsc4py import PETSc
from utils import update_current_density
from utils3D import SupplyCurrentDensity, PMMagnetization, MagneticField3D
from utils2 import DerivedQuantities2D


def AssembleSystem(a, L, bcs, name):
    # if MPI.COMM_WORLD.rank == 0: 
    #     log.log(loglevel, name + ": Assembling LHS Matrix")
    if not bcs:
        A = dolfinx.fem.petsc.assemble_matrix(dolfinx.fem.form(a))    # assemble_matrix(a)
    else:
        A = dolfinx.fem.petsc.assemble_matrix(dolfinx.fem.form(a), bcs) # assemble_matrix(a, bcs)
    A.assemble()

    # if MPI.COMM_WORLD.rank == 0:
    #     log.log(loglevel, name + ": Assembling RHS Vector")
    b = dolfinx.fem.petsc.assemble_vector(dolfinx.fem.form(L))
    dolfinx.fem.petsc.apply_lifting(dolfinx.fem.form(b), [dolfinx.fem.form(a)], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD,
                  mode=PETSc.ScatterMode.REVERSE)
    dolfinx.fem.set_bc(b, bcs)
    return A, b

def SaveSolution(fnc, t, file, name, units):
    # if MPI.COMM_WORLD.rank == 0:
    #     # log.log(loglevel, name + ": Saving Solution")
    file.write_function(fnc, round(t,3))
    
    fncmax, fncmin = fnc.vector.max()[1], fnc.vector.min()[1]
    # if MPI.COMM_WORLD.rank == 0:
    #     log.log(loglevel, name + " Max: " + \
    #         str("{:.3e}".format(fncmax)) + " " + units)
    #     log.log(loglevel, name + " Min: " + \
    #         str("{:.3e}".format(fncmin)) + " " + units)
        

def solve_pmsm(outdir: Path = Path("results"), progress: bool = False, save_output: bool = False):
    """
    Solve the TEAM 30 problem for a single or three phase engine.

    Parameters
    ==========
    outdir
        Directory to put results in

    plot
        Plot torque and voltage over time

    progress
        Show progress bar for solving in time

    save_output
        Save output to bp-files
    """

    # Parameters
    fname = Path("meshes") / "pmesh3D"           # pmesh4_res_0005    # pmsm mesh {pmesh3, pmesh1, pmesh4}
    omega_u: np.float64 = 62.83                     # Angular speed of rotor [rad/s]    # 600 RPM; 1 RPM = 2pi/60 rad/s
    degree: np.int32 = 1                            # Degree of magnetic vector potential functions space (default: 1)
    apply_torque: bool = False                      # Apply external torque to engine (ignore omega) (default: False)
    form_compiler_options: dict = {} 
    jit_parameters: dict = {}

    # Note: model_parameters, domain_parameters and surface_map imported from generate_pmsm_2D script
    # Model parameters for the PMSM model
    mu_0 = 1.25663753e-6
    model_parameters = {
        "mu_0": 1.25663753e-6,      # Relative permability of air [H/m]=[kg m/(s^2 A^2)]
        "freq": 50,                 # Frequency of excitation,
        "J": 1413810.0970277672,     # 3.1e6 * np.sqrt(2),    # [A/m^2] Current density of copper winding
        "mu": {"Cu": 0.999991*mu_0, "Stator": 100*mu_0, "Rotor": 100*mu_0, "Al": 100*mu_0, "Air": mu_0, "AirGap": mu_0, "PM": 1.04457*mu_0},    # "Al": 1.000022*mu_0,    # permability of material
        "sigma": {"Rotor": 2e6, "Al": 2e6, "Stator": 0, "Cu": 0, "Air": 0, "AirGap": 0, "PM": 6.25e5},  # Conductivity 6
        "densities": {"Rotor": 7850, "Al": 7850, "Stator": 0, "Air": 0, "Cu": 0, "AirGap": 0, "PM": 7500}       # [kg/m^3]
    }
    # "Al": 1.000022*mu_0 "Al": 3.72e7 2700
    freq = model_parameters["freq"]             # 50
    T = 0.1 #0.005 0.01 
    dt_ = 0.001 #0.001 0.002 
    omega_J = 2 * np.pi * freq                  # 376.99111843077515

    # Copper wires and PMs are ordered in counter clock-wise order from angle = 0, 2*np.pi/num_segments...
    domains = {"Air": (1,), "AirGap": (2, 3), "Al": (4,), "Rotor": (5, ), 
               "Stator": (6, ), "Cu": (7, 8, 9, 10, 11, 12),
               "PM": (13, 14, 15, 16, 17, 18, 19, 20, 21, 22)}
    
    # Currents mapping to the domain markers of copper
    currents = {7: {"alpha": 1, "beta": 0}, 8: {"alpha": -1, "beta": 2 * np.pi / 3},
                9: {"alpha": 1, "beta": 4 * np.pi / 3}, 10: {"alpha": -1, "beta": 0},
                11: {"alpha": 1, "beta": 2 * np.pi / 3},
                12: {"alpha": -1, "beta": 4 * np.pi / 3}}

    # Marker for facets, and restriction to use in surface integral of airgap
    surface_map: Dict[str, Union[int, str]] = {"Exterior": 1, "MidAir": 2, "restriction": "+"}

    # Read mesh and cell markers
    with io.XDMFFile(MPI.COMM_WORLD, f"{fname}.xdmf", "r") as xdmf:
        mesh = xdmf.read_mesh()                                 
        ct = xdmf.read_meshtags(mesh, name="Cell_markers")      
        tdim = mesh.topology.dim        # 2
        mesh.topology.create_connectivity(tdim - 1, 0)
        ft = xdmf.read_meshtags(mesh, name="Facet_markers")     
    
    # Create DG 0 function for mu and sigma
    DG0 = fem.FunctionSpace(mesh, ("DG", 0))
    mu = fem.Function(DG0)
    sigma = fem.Function(DG0)
    density = fem.Function(DG0)
    for (material, domain) in domains.items():
        for marker in domain:
            cells = ct.find(marker)
            mu.x.array[cells] = model_parameters["mu"][material]        
            sigma.x.array[cells] = model_parameters["sigma"][material]      
            density.x.array[cells] = model_parameters["densities"][material]

    # Define problem function space
    cell = mesh.ufl_cell()                              # cell = 'tetrahedron', degree = 1
    VE = ufl.VectorElement("Lagrange", cell, degree, dim = 3)               
    V_ned = basix.ufl.element("N1curl", mesh.basix_cell(), degree)
    # A_space = fem.functionspace(mesh, nedelec_elem)               
    FE = ufl.FiniteElement("Lagrange", cell, degree)    
    ME = ufl.MixedElement([V_ned, FE])                     
    VQ = fem.FunctionSpace(mesh, ME)        

    # Define test, trial and functions for previous timestep

    Az, V = ufl.TrialFunctions(VQ)
    vz, q = ufl.TestFunctions(VQ)

    # CG_V  = VectorElement("CG", cell, order_v)
    # V_V    = ufl.FunctionSpace(mesh, VE)
    V_V = fem.FunctionSpace(mesh, VE)      # check fem.FunctionSpace(mesh, ("CG", 1, (3,)))
    V_N = fem.functionspace(mesh, V_ned)
    # A0 = Function(V_V)
    AnVn = fem.Function(VQ)
    An, _ = ufl.split(AnVn)  # Solution at previous time step   # may cause error because of 3D function space [Ax, Ay, Az, V]
    
    # J0z = fem.Function(V_V)  # Current density
    J0z = fem.Function(DG0)
    # Supply current (A)
    sp_current  = 35.00
    # Current density magnitude (A/m^2)
    jsource_amp = sp_current/2.47558E-05
    jsexp = SupplyCurrentDensity()
    jsexp.amp = jsource_amp
    jsexp.omega = omega_J
    jsource = fem.Function(V_V)
    jsource.interpolate(jsexp.eval)
    jsource.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT,
                            mode=PETSc.ScatterMode.FORWARD)
        

    # Create integration sets
    Omega_n = domains["Cu"] + domains["Stator"] + domains["Air"] + domains["AirGap"]
    Omega_c = domains["Rotor"] + domains["Al"] + domains["PM"]
    Omega_pm = domains["PM"]
    
    # Magnetization part
    # coercivity = 8.38e5  # [A/m]   
    # DG0v = fem.FunctionSpace(mesh, ("DG", 0, (3,)))     # check
    # Mvec = fem.Function(DG0v)
    # Create magnetization excitation
    # Remanent magnetic flux density (T)
    msource_mag_T = 1.09999682447133
    # Permanent Magnetization (A/m)
    msource_mag   = (msource_mag_T*1e7)/(4*math.pi)
    msexp = PMMagnetization()
    msexp.mag = msource_mag
    msource = fem.Function(V_V)
    msource.interpolate(msexp.eval)
    msource.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT,
                            mode=PETSc.ScatterMode.FORWARD)
    
    # Create integration measures
    dx = ufl.Measure("dx", domain=mesh, subdomain_data=ct)

    # Define temporal and spatial parameters
    dt = fem.Constant(mesh, dt_)
    x = ufl.SpatialCoordinate(mesh)

    omega = fem.Constant(mesh, default_scalar_type(omega_u))

    # Motion voltage term
    (x,y,z) = ufl.SpatialCoordinate(mesh)
    radius  = ufl.as_vector((x,y,0))

    # Magnetization term
    # curl_vz = ufl.as_vector((vz.dx(1), -vz.dx(0)))
    # mag_term =  (mu_0/mu) * ufl.inner( Mvec , curl_vz) * dx(Omega_pm) 
    
    f_a =   + dt / mu * ufl.inner(ufl.grad(Az), ufl.grad(vz)) * dx(Omega_n + Omega_c) \
            + sigma * ufl.inner((Az - An), vz) * dx(Omega_c) \
            + dt * sigma * ufl.inner(ufl.grad(V), vz) * dx(Omega_c) \
            + dt * sigma * ufl.inner(ufl.cross(omega * radius, ufl.curl(Az)), vz) * dx(Omega_c) \
            - dt * J0z * vz[2] * dx(Omega_n) \
            - dt * (mu_0/mu) * ufl.inner( msource , ufl.curl(vz)) * dx(Omega_pm)
            # - dt * ufl.inner(J0z, vz) * dx(Omega_n) \

    f_v =   + dt * sigma * ufl.inner(ufl.grad(V), ufl.grad(q)) * dx(Omega_n + Omega_c) \
            + sigma * ufl.inner((Az - An), ufl.grad(q)) * dx(Omega_c) \
            + dt * sigma * ufl.inner(ufl.cross(omega * radius, ufl.curl(Az)), vz) * dx(Omega_c)

    form_av = f_a + f_v
    a, L = ufl.system(form_av)

    # Find all dofs in Omega_n for Q-space
    cells_n = np.hstack([ct.find(domain) for domain in Omega_n])
    Q, _ = VQ.sub(1).collapse()
    deac_dofs = fem.locate_dofs_topological((VQ.sub(1), Q), tdim, cells_n)

    # Create zero condition for V in Omega_n
    zeroQ = fem.Function(Q)
    zeroQ.x.array[:] = 0
    bc_Q = fem.dirichletbc(zeroQ, deac_dofs, VQ.sub(1))

    # Create external boundary condition for V space
    V_, _ = VQ.sub(0).collapse()
    tdim = mesh.topology.dim
    mesh.topology.create_connectivity(tdim - 1, tdim)
    boundary_facets = dolfinx.mesh.exterior_facet_indices(mesh.topology)
    bndry_dofs = fem.locate_dofs_topological((VQ.sub(0), V_), tdim - 1, boundary_facets)
    zeroV = fem.Function(V_)
    zeroV.x.array[:] = 0
    bc_V = fem.dirichletbc(zeroV, bndry_dofs, VQ.sub(0))
    bcs = [bc_V, bc_Q]

    # Create sparsity pattern and matrix with additional non-zeros on diagonal
    cpp_a = fem.form(a, form_compiler_options=form_compiler_options,
                     jit_options=jit_parameters)
    pattern = fem.create_sparsity_pattern(cpp_a)
    block_size = VQ.dofmap.index_map_bs
    deac_blocks = deac_dofs[0] // block_size
    pattern.insert_diagonal(deac_blocks)
    pattern.finalize()


    # Initialize B
    DG_V  = ufl.VectorElement("DG", cell, 0)
    V_DGV  = dolfinx.fem.FunctionSpace(mesh, DG_V) 
    # Magnetic Flux Density (B)
    u_b, v_b = ufl.TrialFunction(V_DGV), ufl.TestFunction(V_DGV)
    B = dolfinx.fem.Function(V_DGV)
    B.name = "B"
    file_b = dolfinx.io.XDMFFile(mesh.comm, str(outdir / "B.xdmf"), "w")
    file_b.write_mesh(mesh)
    
    # Create matrix based on sparsity pattern
    A = cpp.la.petsc.create_matrix(mesh.comm, pattern)
    A.zeroEntries()
    if not apply_torque:
        A.zeroEntries()
        _petsc.assemble_matrix(A, cpp_a, bcs=bcs)  # type: ignore
        A.assemble()

    # Create inital vector for LHS
    cpp_L = fem.form(L, form_compiler_options=form_compiler_options,
                     jit_options=jit_parameters)
    b = _petsc.create_vector(cpp_L)

    # Create solver
    solver = PETSc.KSP().create(mesh.comm)  # type: ignore
    solver.setOperators(A)
    prefix = "AV_"

    # Give PETSc solver options a unique prefix
    solver_prefix = f"PMSM_solve_{id(solver)}"    
    solver.setOptionsPrefix(solver_prefix)

    # Set PETSc options
    opts = PETSc.Options()  # type: ignore
    opts.prefixPush(solver_prefix)
    # petsc_options: dict = {"ksp_type": "preonly", "pc_type": "lu"}
    # for k, v in petsc_options.items():
    #     opts[k] = v
    opts["ksp_type"] = "gmres"
    # opts["ksp_converged_reason"] = None
    # opts["ksp_monitor_true_residual"] = None
    opts["ksp_type"] = "gmres"
    opts["ksp_gmres_modifiedgramschmidt"] = None
    opts["ksp_diagonal_scale"] = None
    opts["ksp_gmres_restart"] = 500
    opts["ksp_rtol"] = 1e-08
    opts["ksp_max_it"] = 50000
    opts["pc_type"] = "bjacobi"
    # opts["pc_view"] = None
    # opts["ksp_monitor"] = None
    # opts["ksp_view"] = None
    solver.setFromOptions()
    # opts.prefixPop()
    # solver.setFromOptions()
    # solver.setOptionsPrefix(prefix)
    # solver.setFromOptions()

    # # Derived Quantity Solver
    solver_dq = PETSc.KSP().create(mesh.comm)
    solver_dq.setOptionsPrefix("DQ_")
    solver_dq.setType("preonly")
    solver_dq.getPC().setType("lu")
    opts_dq = PETSc.Options("DQ_")
    opts_dq["ksp_type"] = "cg"
    opts_dq["pc_type"] = "bjacobi"
    opts_dq["ksp_monitor"] = None
    opts_dq["ksp_view"] = None
    solver_dq.setFromOptions()



    # Function for containg the solution
    AzV = fem.Function(VQ)
    A_out = fem.Function(V_N)

    A_out = AzV.sub(0).collapse()
    V_out = AzV.sub(1).collapse()

    # Post-processingfunction for projecting the magnetic field potential
    post_B = MagneticField3D(AzV)

    A_out.name = "A"
    post_B.B.name = "B"
    V_out.name = "V"
    
    # Create output file
    if save_output:
        Az_vtx = VTXWriter(mesh.comm, str(outdir / "A.bp"), [A_out])
        B_vtx = VTXWriter(mesh.comm, str(outdir / "B.bp"), [post_B.B])
        V_vtx = VTXWriter(mesh.comm, str(outdir / "V.bp"), [V_out])

    # Computations needed for adding addiitonal torque to engine
    x = ufl.SpatialCoordinate(mesh)
    r = ufl.sqrt(x[0]**2 + x[1]**2)
    L = 1  # Depth of domain
    I_rotor = mesh.comm.allreduce(fem.assemble_scalar(fem.form(L * r**2 * density * dx(Omega_c))))

    # Post proc variables
    num_steps = int(T / float(dt.value))
    torques = np.zeros(num_steps + 1, dtype=default_scalar_type)
    torques_vol = np.zeros(num_steps + 1, dtype=default_scalar_type)
    times = np.zeros(num_steps + 1, dtype=default_scalar_type)
    omegas = np.zeros(num_steps + 1, dtype=default_scalar_type)
    omegas[0] = omega_u
    pec_tot = np.zeros(num_steps + 1, dtype=default_scalar_type)
    pec_steel = np.zeros(num_steps + 1, dtype=default_scalar_type)
    VA = np.zeros(num_steps + 1, dtype=default_scalar_type)
    VmA = np.zeros(num_steps + 1, dtype=default_scalar_type)
    # Generate initial electric current in copper windings
    t = 0.
    update_current_density(J0z, omega_J, t, ct, currents)
    # update_magnetization(Mvec, coercivity, omega_u, t, ct, domains, pm_orientation)
    if MPI.COMM_WORLD.rank == 0 and progress:
        progressbar = tqdm.tqdm(desc="Solving time-dependent problem",
                                total=int(T / float(dt.value)))

    for i in range(num_steps):
        # Update time step and current density
        if MPI.COMM_WORLD.rank == 0 and progress:
            progressbar.update(1)
        t += float(dt.value)
        jsexp.t = t
        jsource.interpolate(jsexp.eval)
        jsource.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT,
                                mode=PETSc.ScatterMode.FORWARD)
        update_current_density(J0z, omega_J, t, ct, currents)
        # update_magnetization(Mvec, coercivity, omega_u, t, ct, domains, pm_orientation)

        # Reassemble RHS
        with b.localForm() as loc_b:
            loc_b.set(0)
        _petsc.assemble_vector(b, cpp_L)
        _petsc.apply_lifting(b, [cpp_a], [bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)  # type: ignore
        fem.set_bc(b, bcs)

        # Solve problem
        solver.solve(b, AzV.vector)
        AzV.x.scatter_forward()

        # Update previous time step
        AnVn.x.array[:] = AzV.x.array
        AnVn.x.scatter_forward()
        
        # Update B
        A_res = AzV.sub(0).collapse()
        form_b = + ufl.inner(u_b,v_b)*dx \
                - ufl.inner(ufl.curl(A_res),v_b)*dx \

        a_b, L_b = ufl.system(form_b)

        A_b, b_b = AssembleSystem(a_b, L_b, [], "B")
        
        A_b = dolfinx.fem.petsc.assemble_matrix(dolfinx.fem.form(a_b)) #assemble_matrix(a_b)
        A_b.assemble()
        
        # if MPI.COMM_WORLD.rank == 0:
        #     log.log(loglevel, "B: Calculating")
        solver_dq.setOperators(A_b)
        solver_dq.solve(b_b, B.vector)
        # t_08.stop()

        # t_09 = Timer("09 Save B Solution")
        SaveSolution(B, t, file_b, "B", "Tm")


        # Update rotational speed 
        # omegas[i + 1] = float(omega.value)
        print("Az = ", min(AzV.sub(0).collapse().x.array[:]), max(AzV.sub(0).collapse().x.array[:]))
        print("B_val = ", min(post_B.B.x.array[:]), max(post_B.B.x.array[:]))
        # Write solution to file
        if save_output:
            post_B.interpolate()
            A_out.x.array[:] = AzV.sub(0).collapse().x.array[:]
            V_out.x.array[:] = AzV.sub(1).collapse().x.array[:]
            Az_vtx.write(t)
            B_vtx.write(t)
            V_vtx.write(t)
    b.destroy()
    # print("B_val = ", post_B.B.x.array[:])
    if save_output:
        Az_vtx.close()
        B_vtx.close()
        V_vtx.close()

    elements = mesh.topology.index_map(mesh.topology.dim).size_global
    num_dofs = VQ.dofmap.index_map.size_global * VQ.dofmap.index_map_bs
    # Print values for last period
    if mesh.comm.rank == 0:
        print(f"{omega_u}, {freq}, {degree}, {elements}, {num_dofs} ")


if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser(
        description="Script to  solve the PMSM problem",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--progress', dest='progress', action='store_true',
                        help="Show progress bar", default=False)
    parser.add_argument('--output', dest='output', action='store_true',
                        help="Save output to VTXFiles files", default=False)

    args = parser.parse_args()
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime("%b_%d_%H_%M_%S")
    outdir = Path(f"PMSM_3D{formatted_datetime}")
    outdir.mkdir(exist_ok=True)
    print(f"Saving to PMSM_3D{formatted_datetime}")
    solve_pmsm(outdir=outdir, progress=args.progress, save_output=args.output)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")