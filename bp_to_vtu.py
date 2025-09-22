import adios2
import vtk
import numpy as np
from vtk.util import numpy_support
import argparse
import os

class BPToVTUConverter:
    def __init__(self):
        self.reader = None
        self.writer = None
    
    def read_bp_file(self, bp_filename):
        """
        Read BP file using ADIOS2
        """
        try:
            # Initialize ADIOS2
            adios = adios2.ADIOS()
            io = adios.DeclareIO("reader")
            
            # Open BP file for reading
            engine = io.Open(bp_filename, adios2.Mode.Read)
            
            # Get available variables
            variables = io.AvailableVariables()
            print(f"Available variables in {bp_filename}:")
            for name, info in variables.items():
                print(f"  {name}: {info}")
            
            data = {}
            
            # Read all variables
            for var_name in variables:
                var = io.InquireVariable(var_name)
                if var:
                    # Get variable shape and type
                    shape = var.Shape()
                    var_type = var.Type()
                    print(f"Reading {var_name} with shape {shape} and type {var_type}")
                    
                    # Read the variable data
                    engine.BeginStep()
                    engine.Get(var, adios2.Mode.Sync)
                    engine.EndStep()
                    
                    # Store data
                    data[var_name] = np.array(var)
            
            engine.Close()
            return data
            
        except Exception as e:
            print(f"Error reading BP file: {e}")
            return None
    
    def create_vtu_from_data(self, data, output_filename):
        """
        Create VTU file from data dictionary
        """
        try:
            # Create unstructured grid
            ugrid = vtk.vtkUnstructuredGrid()
            
            # Assume we have coordinates (modify based on your data structure)
            coords_keys = ['x', 'y', 'z', 'coordinates', 'points']
            coord_data = None
            
            for key in coords_keys:
                if key in data:
                    coord_data = data[key]
                    break
            
            if coord_data is None:
                # Try to find coordinate-like data
                for key, value in data.items():
                    if 'coord' in key.lower() or 'point' in key.lower():
                        coord_data = value
                        break
            
            if coord_data is not None:
                # Handle different coordinate formats
                if coord_data.ndim == 1:
                    # Assume 1D array needs to be reshaped
                    n_points = len(coord_data) // 3
                    if len(coord_data) % 3 == 0:
                        coord_data = coord_data.reshape(n_points, 3)
                
                if coord_data.ndim == 2 and coord_data.shape[1] == 3:
                    # Create VTK points
                    points = vtk.vtkPoints()
                    for i in range(coord_data.shape[0]):
                        points.InsertNextPoint(coord_data[i, 0], coord_data[i, 1], coord_data[i, 2])
                    
                    ugrid.SetPoints(points)
                    
                    # Create cells (assuming point cloud, you may need to modify for your data)
                    for i in range(points.GetNumberOfPoints()):
                        ugrid.InsertNextCell(vtk.VTK_VERTEX, 1, [i])
            
            # Add scalar data as point data
            point_data = ugrid.GetPointData()
            
            for var_name, var_data in data.items():
                if var_name not in coords_keys and 'coord' not in var_name.lower():
                    # Convert numpy array to VTK array
                    if var_data.ndim == 1:
                        vtk_array = numpy_support.numpy_to_vtk(var_data)
                        vtk_array.SetName(var_name)
                        point_data.AddArray(vtk_array)
                        
                        # Set first scalar array as active
                        if point_data.GetNumberOfArrays() == 1:
                            point_data.SetActiveScalars(var_name)
            
            # Write VTU file
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(output_filename)
            writer.SetInputData(ugrid)
            writer.SetDataModeToAscii()  # Use SetDataModeToBinary() for binary output
            writer.Write()
            
            print(f"Successfully wrote VTU file: {output_filename}")
            return True
            
        except Exception as e:
            print(f"Error creating VTU file: {e}")
            return False
    
    def convert(self, bp_filename, vtu_filename):
        """
        Main conversion function
        """
        print(f"Converting {bp_filename} to {vtu_filename}")
        
        # Read BP file
        data = self.read_bp_file(bp_filename)
        if data is None:
            return False
        
        # Create VTU file
        return self.create_vtu_from_data(data, vtu_filename)

def main():
    parser = argparse.ArgumentParser(description='Convert BP file to VTU format')
    parser.add_argument('input', help='Input BP file path')
    parser.add_argument('output', help='Output VTU file path')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found")
        return
    
    # Create converter and run conversion
    converter = BPToVTUConverter()
    success = converter.convert(args.input, args.output)
    
    if success:
        print("Conversion completed successfully!")
    else:
        print("Conversion failed!")

# Example usage as a script
if __name__ == "__main__":
    # If no command line arguments, show example usage
    import sys
    if len(sys.argv) == 1:
        print("Example usage:")
        print("python bp_to_vtu.py input.bp output.vtu")
        print("\nOr use as a module:")
        print("converter = BPToVTUConverter()")
        print("converter.convert('input.bp', 'output.vtu')")
    else:
        main()

# Alternative simple function for direct use
def convert_bp_to_vtu(bp_file, vtu_file):
    """
    Simple function to convert BP to VTU
    
    Args:
        bp_file (str): Path to input BP file
        vtu_file (str): Path to output VTU file
    
    Returns:
        bool: True if successful, False otherwise
    """
    converter = BPToVTUConverter()
    return converter.convert(bp_file, vtu_file)