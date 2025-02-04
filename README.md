# PMSM-FEniCS

This repository contains FEniCS implementation of a Permanent Magnet Synchronous Motor (PMSM) Model in 2D and 3D. [(Reference)](https://doi.org/10.1016/j.finel.2022.103755)

## 2D modelling

- `generate_pmsm_2D.py`: A script that generates the PMSM 2D meshes and saves them to xdmf format. 
- `pmsm.py`: Script for solving the PMSM 2D model.  
- `utils.py`: File containing utillity functions used in pmsm.py
- `pmsm_msh.py`: Script for running different mesh resolutions. 
- `pmsm_tor.py`: Script for calculating torques. 
- `2D results`: Folder containing 2D Model results
- `meshes`: Folder containing different meshes used for simulation
- `mesh images`: Folder containing images of meshes used for simulation

## 3D modelling

- `generate_pmsm_3D.py`: A script that generates the PMSM 3D meshes and saves them to xdmf format. 
- `pmsm_3D.py`: Script for solving the PMSM 2D model.
- `utils3D.py`: File containing utillity functions used in pmsm_3D.py
- `3D results`: Folder containing 3D Model results
---
## 2D Model Results
- [Gmsh](https://gmsh.info/) - Used for viewing and verifying the mesh.
- [ParaView](https://www.paraview.org/) - Used for analyzing and visualizing the results.
<table>
  <tr>
    <td align="center">
      <a href="./meshes/mesh%20screenshots/2D%20mesh%20-%20res%200.001.png">
        <img src="./meshes/mesh%20screenshots/2D%20mesh%20-%20res%200.001.png" alt="PMSM 2D Mesh" width="300">
      </a>
      <p><b>PMSM 2D Mesh</b></p>
    </td>
    <td align="center">
      <a href="./2D%20results/PMSM2D_Az.png">
        <img src="./2D%20results/PMSM2D_Az.png" alt="Magnetic Vector Potential" width="300">
      </a>
      <p><b>Magnetic Vector Potential</b></p>
      <p><a href="https://drive.google.com/file/d/1GrcKroc-dno4-W_8fjqnpYRnWG8lCafS/view?usp=sharing">Watch it rotate</a></p>
    </td>
    <td align="center">
      <a href="./2D%20results/PMSM2D_B.png">
        <img src="./2D%20results/PMSM2D_B.png" alt="Magnetic Flux Density" width="300">
      </a>
      <p><b>Magnetic Flux Density</b></p>
      <p><a href="https://drive.google.com/file/d/1d07aff6dNZa5njJSuvVtAdGsNQwggtUH/view?usp=sharing">Watch it rotate</a></p>
    </td>
  </tr>
</table>

## 3D Model Results
### 3D Mesh
<table>
  <tr>
    <td align="center">
      <a href="./meshes/mesh%20screenshots/3D_mesh_1.png">
        <img src="./meshes/mesh%20screenshots/3D_mesh_1.png" alt="PMSM 3D Mesh" width="300">
      </a>
      <p><b>Isometric View</b></p>
    </td>
    <td align="center">
      <a href="./meshes/mesh%20screenshots/3D_mesh_2.png">
        <img src="./meshes/mesh%20screenshots/3D_mesh_2.png" alt="PMSM 3D Mesh" width="300">
      </a>
      <p><b>Side View</b></p>
    </td>
    <td align="center">
      <a href="./meshes/mesh%20screenshots/3D_mesh_3.png">
        <img src="./meshes/mesh%20screenshots/3D_mesh_3.png" alt="PMSM 3D Mesh" width="300">
      </a>
      <p><b>Top View</b></p>
    </td>
  </tr>
</table>

***
## Installation
This repository is based on the `TEAM30` project and is configured to run inside a Docker container with `dolfinx`. This guide provides step-by-step instructions to set up and run the project.


### **Prerequisites**

Ensure you have the following installed on your system:

- [Docker](https://docs.docker.com/get-docker/) (Required for running the containerized environment)
- [Git](https://git-scm.com/downloads) (Required for cloning the repository)
- [VS Code](https://code.visualstudio.com/download) (Recommended for editing and debugging the project)

### **Setup**

#### **Step 1: Clone the Repository**
```bash
git clone git@github.com:abhinavtk7/pmsm-fenics.git
cd pmsm-fenics
```

#### **Step 2: Start the Docker Container**
Run the following command from the project root to start the container:
```bash
docker run -ti -v $(pwd):/root/shared -w /root/shared/ --shm-size=512m --name=pmsm ghcr.io/fenics/dolfinx/dolfinx:v0.7.0
```

#### **Step 3: Install Dependencies**
Inside the container, install required Python packages:
```bash
python3 -m pip install tqdm pandas
```

#### **Step 4: Verify Installation**
Run the python code `pmsm_msh.py` with resolution of 0.01 mm inside the container to have a quick check that everything is working correctly:
```bash
python3 pmsm_msh.py --res 3 --progress
```
---
### **Usage**
#### **Starting the Docker Container Again**
If you exit the container, restart it using:
```bash
docker start -ai pmsm
```
If the container was removed, recreate it using the command from Step 2.
#### **Stopping the Docker Container**
Exit the container by running:
```bash
exit
```
Or stop it from another terminal:
```bash
docker stop pmsm
```
#### **Removing the Docker Container**
To remove the container completely:
```bash
docker rm pmsm
```
### **Using VS Code for editing and debugging**
If you want to edit or debug the code, you can do it using VS Code.
##### **Open VS Code and Attach to the Running Container**
- Launch VS Code.
- If you haven't installed "Remote - Containers" Extension, install it.
  - Go to Extensions (Ctrl + Shift + X)
  - Search for "Remote - Containers" and click Install. 
- Press Ctrl + Shift + P to open the command palette.
- Type and select "Attach to Running Container".
- Choose the container named `pmsm` (Make sure you've started the container).
- Once attached, you will be inside the `pmsm` container.
