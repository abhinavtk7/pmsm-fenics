# Reference code: https://github.com/Wells-Group/TEAM30

from typing import Dict, Tuple
import basix.ufl
import numpy as np
import ufl
from dolfinx import cpp, default_scalar_type, fem
from mpi4py import MPI
from petsc4py import PETSc
import math
from scipy.spatial.transform import Rotation

domains = {"Air": (1,), "AirGap": (2, 3), "Al": (4,), "Rotor": (5, ), 
               "Stator": (6, ), "Cu": (7, 8, 9, 10, 11, 12),
               "PM": (13, 14, 15, 16, 17, 18, 19, 20, 21, 22)}

currents = {7: {"alpha": 1, "beta": 0}, 8: {"alpha": -1, "beta": 2 * np.pi / 3},
                9: {"alpha": 1, "beta": 4 * np.pi / 3}, 10: {"alpha": -1, "beta": 0},
                11: {"alpha": 1, "beta": 2 * np.pi / 3},
                12: {"alpha": -1, "beta": 4 * np.pi / 3}}


from dolfinx.geometry import bb_tree
from dolfinx.geometry import compute_collisions_points, compute_colliding_cells
from mpi4py import MPI
from dolfinx import cpp, fem, io, default_scalar_type
from pathlib import Path


fname = Path("meshes") / "pmesh3D" 
with io.XDMFFile(MPI.COMM_WORLD, f"{fname}.xdmf", "r") as xdmf:
        mesh = xdmf.read_mesh()                                 
        ct = xdmf.read_meshtags(mesh, name="Cell_markers")      
        tdim = mesh.topology.dim        # 3
        mesh.topology.create_connectivity(tdim - 1, 0)
        ft = xdmf.read_meshtags(mesh, name="Facet_markers") 

# Build a bounding box tree for the mesh
tree = bb_tree(mesh, mesh.topology.dim)

# Function to find the cell containing a point (x, y, z)
def find_cell_from_point(x, y, z):
    point = np.array([x, y, z], dtype=np.float64)
    
    # Query the bounding box tree
    cell_candidates = compute_collisions_points(tree, point)
    colliding_cells = compute_colliding_cells(mesh, cell_candidates, point)
    
    # Check if any cell was found
    if len(colliding_cells) > 0:
        return colliding_cells[0]  # Return the first cell that contains the point
    else:
        return None

class SupplyCurrentDensity:
    def __init__(self):
        self.cw = 1.
        self.amp = 1.
        self.omega = 1.
        self.t = 0.
        
        self.yl = 3.2223088765338/1000.
        self.zl = 40./1000.
        self.m = 0.577

    def eval(self, coord):
        """PMSM Current Density Excitation"""
        values = np.zeros((3, coord.shape[1]))
        
        yl, zl, m = self.yl, self.zl, self.m
        cw, amp, omega, t = self.cw, self.amp, self.omega, self.t
        
        phase_a = amp*math.sin( omega*t )
        phase_b = amp*math.sin( omega*t + (2*math.pi/3) )
        phase_c = amp*math.sin( omega*t - (2*math.pi/3) )
        

        for i in range(coord.shape[1]):
            x, y, z = coord[0][i], coord[1][i], coord[2][i] 
            cell_number = find_cell_from_point(x, y, z)
            marker = ct.values[np.where(ct.indices == cell_number)[0][0]]
            if marker not in [7, 8, 9, 10, 11, 12]:
                continue
            # print(cell_number, marker)
            if marker == 7:     # S0
                rotangle = 0
                pmag = phase_a
                coilcw = -1
            elif marker == 8:   # S1
                rotangle = 60
                pmag = phase_b
                coilcw = -1
            elif marker == 9:   # S2
                rotangle = 120
                pmag = phase_c
                coilcw = 1
            elif marker == 10:     # S3
                rotangle = 180
                pmag = phase_a
                coilcw = 1
            elif marker == 11:   # S4
                rotangle = 240
                pmag = phase_b
                coilcw = 1
            elif marker == 12:   # S5
                rotangle = 300
                pmag = phase_c
                coilcw = -1
            
            # rotate coordinates to S0
            r = Rotation.from_euler('z', -1*rotangle, degrees=True)
            vec = r.apply([x, y, z])
            x, y, z = vec[0], vec[1], vec[2]
            
            # calculate which region in S0 the point is in
            # R1 (linear)
            if (y <= yl and y >= -1*yl) and (z > zl):
                values[1, i] = cw*coilcw*-1*pmag
                 
            # R2 (rotation)
            elif y < -1*yl and z > zl:
                dy, dz = y + yl, z - zl
                r = math.sqrt(dy**2 + dz**2)
                values[1, i] = cw*coilcw*-1*pmag*(dz/r)
                values[2, i] = cw*coilcw*pmag*(dy/r)
                
            # R3 (linear)
            elif (z <= zl and z >= -1*zl) and (y < -1*yl):
                values[2, i] = cw*coilcw*-1*pmag
                
            # R4 (rotation)
            elif y < -1*yl and z < -1*zl:
                dy, dz = y + yl, z + zl
                r = math.sqrt(dy**2 + dz**2)
                values[1, i] = cw*coilcw*-1*pmag*(dz/r)
                values[2, i] = cw*coilcw*pmag*(dy/r)
                
            # R5 (linear)
            elif (y <= yl and y >= -1*yl) and (z < -1*zl):
                values[1, i] = cw*coilcw*pmag
                
            # R6 (rotation)
            elif y > yl and z < -1*zl:
                dy, dz = y - yl, z + zl
                r = math.sqrt(dy**2 + dz**2)
                values[1, i] = cw*coilcw*-1*pmag*(dz/r)
                values[2, i] = cw*coilcw*pmag*(dy/r)
                
            # R7 (linear)
            elif (z <= zl and z >= -1*zl) and (y > yl):
                values[2, i] = cw*coilcw*pmag
                
            # R8 (rotation)
            elif y > yl and z > zl:
                dy, dz = y - yl, z - zl
                r = math.sqrt(dy**2 + dz**2)
                values[1, i] = cw*coilcw*-1*pmag*(dz/r)
                values[2, i] = cw*coilcw*pmag*(dy/r)
            
            # rotate the excitation back to original sector
            rot = Rotation.from_euler('z', rotangle, degrees=True)
            
            vec = [values[0, i], values[1, i], values[2, i]]
            values[0, i] = rot.apply(vec)[0]
            values[1, i] = rot.apply(vec)[1]
            values[2, i] = rot.apply(vec)[2]
            
        return values


def update_magnetization(Mvec, coercivity, omega_u, t, ct, domains, pm_orientation):
        block_size = 3 # Mvec.function_space.dofmap.index_map_bs = 2 for 2D; 3 for 3D
        coercivity = 8.38e5  # [A/m]   
        sign = 1

        for (material, domain) in domains.items():
            if material == 'PM':
                for marker in domain:
                    inout = 1 if marker in [13, 15, 17, 19, 21] else -1
                    angle = pm_orientation[marker] + omega_u * t
                    Mx = coercivity * np.cos(angle) * sign * inout
                    My = coercivity * np.sin(angle) * sign * inout
                    Mz = 0.0
                    cells = ct.find(marker)
                    for cell in cells:
                        idx = block_size * cell
                        Mvec.x.array[idx + 0] = Mx
                        Mvec.x.array[idx + 1] = My
                        Mvec.x.array[idx + 2] = Mz

        # Mvec.x.scatter_forward()
        Mvec.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT,
                                mode=PETSc.ScatterMode.FORWARD)



class PMMagnetization:
    def __init__(self):
        self.sign = 1.
        self.mag = 1.
        
    def eval(self, coord):
        """PMSM Magnetization Excitation"""
        values = np.zeros((3, coord.shape[1]))
        
        mag, sign = self.mag, self.sign
        pm_spacing = (np.pi / 6) + (np.pi / 30)
        pm_angles = np.asarray([i * pm_spacing for i in range(10)], dtype=np.float64)
        # [13, 14, 15, 16, 21, 22]
        # apply load for aligned domain
        for i in range(coord.shape[1]):
            x, y, z = coord[0][i], coord[1][i], coord[2][i] 
            cell_number = find_cell_from_point(x, y, z)
            marker = ct.values[np.where(ct.indices == cell_number)[0][0]]
            if marker not in [ 17, 18, 19, 20]:
                continue
            # print(cell_number, marker)
            if marker == 13: 
                rotangle = pm_angles[0]
                inout = 1
            elif marker == 14: 
                rotangle = pm_angles[1]
                inout = -1
            elif marker == 15: 
                rotangle = pm_angles[2]
                inout = 1
            elif marker == 16: 
                rotangle = pm_angles[3]
                inout = -1
            elif marker == 17: 
                rotangle = pm_angles[4]
                inout = 1
            elif marker == 18: 
                rotangle = pm_angles[5]
                inout = -1
            elif marker == 19: 
                rotangle = pm_angles[6]
                inout = 1
            elif marker == 20: 
                rotangle = pm_angles[7]
                inout = -1
            elif marker == 21: 
                rotangle = pm_angles[8]
                inout = 1
            elif marker == 22: 
                rotangle = pm_angles[9]
                inout = -1

            # create unit vector in correct direction
            r = Rotation.from_euler('z', rotangle, degrees=True)
            unitvec = r.apply([1, 0, 0])
            result = sign*mag*inout*unitvec
            
            values[0, i] = result[0]
            values[1, i] = result[1]
            values[2, i] = result[2]
            
        return values
    

class MagneticField3D():
    def __init__(self, AzV: fem.Function,
                 form_compiler_options: dict = {}, jit_parameters: dict = {}):
        """
        Class for interpolate the magnetic vector potential (here as the first part of the mixed function AvZ)
        to the magnetic flux intensity B=curl(A)

        Parameters
        ==========
        AzV
            The mixed function of the magnetic vector potential Az and the Scalar electric potential V
        """
        degree = AzV.function_space.ufl_element().degree()
        mesh = AzV.function_space.mesh
        cell = mesh.ufl_cell()

        # Create dolfinx Expression for electromagnetic field B (post processing)
        # Use minimum DG 1 as VTXFile only supports CG/DG>=1
        # el_B = basix.ufl.element("DG", cell.cellname(),
        #                          max(degree - 1, 1),
        #                          shape=(mesh.geometry.dim,),
        #                          gdim=mesh.geometry.dim)
        VE = ufl.VectorElement("Lagrange", cell, degree, dim = 3) 
        VB = fem.FunctionSpace(mesh, VE)  # el_B
        self.B = fem.Function(VB)

        A_vector = AzV.sub(0).collapse()
        B_3D = ufl.curl(A_vector)     # ufl.as_vector((AzV[0].dx(1), -AzV[0].dx(0)))
        self.Bexpr = fem.Expression(B_3D, VB.element.interpolation_points(),
                                    form_compiler_options=form_compiler_options,
                                    jit_options=jit_parameters)

    def interpolate(self):
        """
        Interpolate magnetic field
        """
        self.B.interpolate(self.Bexpr)
