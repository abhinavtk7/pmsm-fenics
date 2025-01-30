# PMSM-FEniCS

This repository contains a DOLFINx implementation of the [PMSM model]([http://www.compumag.org/jsite/images/stories/TEAM/problem30a.pdf](https://doi.org/10.1016/j.finel.2022.103755)).

## 2D modelling

- `generate_pmsm_2D.py`: A script that generates the PMSM 2D meshes and saves them to xdmf format. 
- `pmsm.py`: Script for solving the PMSM 2D model 

## 3D modelling

- `scriptname.py`: description
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

## **Usage**

### **Starting the Docker Container Again**
If you exit the container, restart it using:
```bash
docker start -ai pmsm
```
If the container was removed, recreate it using the command from **Step 2**.

### **Stopping the Docker Container**
Exit the container by running:
```bash
exit
```
Or stop it from another terminal:
```bash
docker stop pmsm
```

### **Removing the Docker Container**
To remove the container completely:
```bash
docker rm pmsm
```

