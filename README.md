# PMSM-FEniCS

This repository contains a DOLFINx implementation of the [PMSM model]([http://www.compumag.org/jsite/images/stories/TEAM/problem30a.pdf](https://doi.org/10.1016/j.finel.2022.103755)).

## 2D modelling

- `generate_pmsm_2D.py`: A script that generates the PMSM 2D meshes and saves them to xdmf format. 
- `pmsm.py`: Script for solving the PMSM 2D model.  
- `utils.py`: File containing utillity functions used in pmsm.py
- `pmsm_msh.py`: 
- `pmsm_tor.py`: 
- `2D results`: Folder containing 2D Model results
- `meshes`: Folder containing different meshes used for simulation
- `mesh images`: Folder containing images of meshes used for simulation

## 3D modelling

- `generate_pmsm_3D.py`: A script that generates the PMSM 3D meshes and saves them to xdmf format. 
- `pmsm_3D.py`: Script for solving the PMSM 2D model.
- `utils3D.py`: File containing utillity functions used in pmsm_3D.py
- `3D results`: Folder containing 3D Model results

## 2D Model Results
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
      <p><a href="https://drive.google.com/file/d/1GrcKroc-dno4-W_8fjqnpYRnWG8lCafS/view?usp=sharing">Watch 'em rotate!</a></p>
    </td>
    <td align="center">
      <a href="./2D%20results/PMSM2D_B.png">
        <img src="./2D%20results/PMSM2D_B.png" alt="Magnetic Flux Density" width="300">
      </a>
      <p><b>Magnetic Flux Density</b></p>
      <p><a href="https://drive.google.com/file/d/1d07aff6dNZa5njJSuvVtAdGsNQwggtUH/view?usp=sharing">Watch 'em rotate!</a></p>
    </td>
  </tr>
</table>


## Installation
This repository is based on the `TEAM30` project and is configured to run inside a Docker container with `dolfinx`. This guide provides step-by-step instructions to set up and run the project.


### **Prerequisites**

Ensure you have the following installed on your system:

- [Docker](https://docs.docker.com/get-docker/) (Required for running the containerized environment)
- [Git](https://git-scm.com/downloads) (Required for cloning the repository)

### **Setup**

#### **Step 1: Clone the Repository**
```bash
git clone git@github.com:abhinavtk7/pmsm-fenics.git
cd pmsm-fenics
```

### **Step 2: Start the Docker Container**
Run the following command from the project root to start the container:
```bash
docker run -ti -v $(pwd):/root/shared -w /root/shared/ --shm-size=512m --name=pmsm ghcr.io/fenics/dolfinx/dolfinx:v0.7.0
```

This command:
- Mounts the current directory to `/root/shared/` inside the container.
- Sets the working directory to `/root/shared/`.
- Allocates shared memory (`--shm-size=512m`) to prevent memory issues.
- Names the container `pmsm`.

### **Step 3: Install Dependencies**
Inside the container, install required Python packages:
```bash
python3 -m pip install tqdm pandas
```

### **Step 4: Verify Installation**
Run the following command inside the container to check that everything is working correctly:
```bash
python3 -m pytest -xvs .
```

If all tests pass, the setup is successful.

### **Removing the Docker Container**
To remove the container completely:
```bash
docker rm pmsm
```
