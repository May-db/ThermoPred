from morfeus import read_xyz, Density
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import streamlit as st
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors
from stmol import showmol
import py3Dmol
from pathlib import Path
import pandas as pd
import os
from streamlit_ketcher import st_ketcher
from rdkit.Chem import rdFingerprintGenerator
import mols2grid
import streamlit.components.v1 as components
import plotly.figure_factory as ff
from typing import Tuple, List

def generate_3D(smiles):
    "Generate 3D coordinates from smiles"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol)
    molstring= Chem.MolToMolBlock(mol)
    return molstring

def visualize_3D(molstring, title="Molecule"):
    "Visualize the molecule in 3D using stmol"
    w, h = 400, 400
    xyzview = py3Dmol.view(width=w,height=w)
    xyzview.addModel(molstring,'mol')
    xyzview.setStyle({'sphere':{'colorscheme':'cyanCarbon', 'scale':0.25}, 'stick':{'colorscheme':'cyanCarbon'}})
    xyzview.zoomTo()
    xyzview.spin()
    xyzview.setBackgroundColor('white')
    with st.container():
        st.markdown(f"### {title} :")
        showmol(xyzview, height = w,width=w)






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

st.title('Create your reaction')
st.caption("Practical Proramming In Chemistry Project")
st.markdown("Draw two molecules and see if the product is stable enough")


# --- Gestion du reset ---
if "mol1" not in st.session_state:
    st.session_state.mol1 = ""
if "mol2" not in st.session_state:
    st.session_state.mol2 = ""


# --- Draw Molecule 1 ---
st.subheader("Draw Molecule 1")
mol1 = st_ketcher(st.session_state.mol1, key="mol1_ketcher", height=400)
if mol1:
    st.session_state.mol1 = mol1
    with st.expander("SMILES Molecule 1"):
        st.code(mol1)
    mol1_3D = generate_3D(mol1)
    if mol1_3D:
        with st.expander("Visualisation 3D Molécule 1", expanded=False):
            visualize_3D(mol1_3D, title="Molécule 1")
    
# --- Draw Molecule 2 ---
with st.container():
    st.subheader("Draw Molecule 2")

    if "show_mol2" not in st.session_state:
        st.session_state.show_mol2 = False

    if st.button("Start to draw molecule 2"):
        st.session_state.show_mol2 = True

    if st.session_state.show_mol2:
        mol2 = st_ketcher(st.session_state.mol2, key="mol2_ketcher", height=400)
        if mol2:
            st.session_state.mol2 = mol2
            with st.expander("SMILES Molecule 2"):
                st.code(mol2)
            mol2_3D = generate_3D(mol2)
            if mol2_3D:
                with st.expander("Visualisation 3D Molecule 2", expanded=False):
                    visualize_3D(mol2_3D, title="Molecule 2")