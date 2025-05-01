from morfeus import read_xyz, Density
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def GeomOptxyz_Energy (molecule_path_xyz: str):
    #test that the file containing the 3D molecule does exist
    if not os.path.exists(molecule_path_xyz):
        raise FileNotFoundError(f"{molecule_path_xyz} not found")
    #does the geometry optimization
    xtb_cmd = ["xtb", molecule_path_xyz, "--opt"]
    #attributes the outputs of the previous command
    output_xyz = "xtbopt.xyz"
    output_log = "xtb.out"
    atoms, coords = read_xyz("xtbopt.xyz")
    result = subprocess.run(xtb_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    #test that the optimization did occur
    if result.returncode != 0:
        print("xTB failed:")
        print(result.stderr)
        raise RuntimeError("xTB optimization failed")
    #gets the energy in the output file (as it is automatically calculated with the optimization)
    energy_hartree = None
    with open(output_log, "r") as f:
        for line in f:
            if "TOTAL ENERGY" in line:
                parts = line.strip().split()
                energy_hartree = float(parts[-2])
                break
    
    return energy_hartree

def Energy_comparison (E_reagent1: float, E_reagent2: float, E_Prod:float):
    #define the energy of the starting point and the variation of energy
    E_start= E_reagent1 + E_reagent2
    dE= E_start - E_Prod
    #setting of the thermodynamic conditions
    if dE>0:
        print ("The reaction is thermodinamically stable at 0 K")
        return True
    if dE==0:
        print("The reaction leads to a thermodinamic equilibrium at 0 K")
    else: 
        print("the reaction is not thermodinamically stable at 0 K")
        return False
    
def smiles_to_3d(smiles, add_H=True, optimize=True):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    if add_H:
        mol=Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())  
    if optimize:
        AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coordinates = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
            
    return elements, coordinates

def write_xyz_file(elements, coordinates, filename):
    with open(filename, 'w') as f:
        f.write(f"{len(elements)}\n\n")
        for element, coord in zip(elements, coordinates):
            f.write(f"{element}{coord[0]:.6f}{coord[1]:.6f} {coord[2]:.6f}\n")