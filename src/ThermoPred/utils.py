from morfeus import read_xyz
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import streamlit as st
from stmol import showmol
import py3Dmol
from streamlit_ketcher import st_ketcher


def generate_3D(smiles):
    "Generate 3D coordinates from smiles"
    mol = Chem.MolFromSmiles(smiles) #convert smiles to an RDKit molecule object
    if mol is None:
        return None
    mol = Chem.AddHs(mol) # add explicit hydrogens
    params = AllChem.ETKDGv3() #use ETKDGv3 to generate 3D coordinates
    params.randomSeed = 42 #reproducibility
    AllChem.EmbedMolecule(mol) #generate 3D coordonates
    molstring= Chem.MolToMolBlock(mol) #convert the 3D molecule to MolBlock format
    return molstring


def visualize_3D(molstring):
    "Visualize the molecule in 3D using stmol and py3Dmol"
    w, h = 400, 400  # Set width and height for the 3D viewer
    xyzview = py3Dmol.view(width=w, height=h)  # Create a 3Dmol.js viewer object
    xyzview.addModel(molstring, 'mol')  # Load the MolBlock string into the viewer as a molecule
    # Set visual style: use spheres and sticks with a cyan color scheme
    xyzview.setStyle({
        'sphere': {'colorscheme': 'cyanCarbon', 'scale': 0.25},
        'stick': {'colorscheme': 'cyanCarbon'}
    })
    xyzview.zoomTo()  # Adjust the camera to fit the molecule in view
    xyzview.spin()  # Enable spinning animation
    xyzview.setBackgroundColor('white')  # Set the background color of the viewer
    with st.container():
        showmol(xyzview, height=w, width=w)  # Render the 3D molecule in Streamlit




def GeomOptxyz_Energy (molecule_path_xyz: str):
    #test that the file containing the 3D molecule does exist

    if not os.path.exists(molecule_path_xyz):
        raise FileNotFoundError(f"Input file {molecule_path_xyz} not found")
    #does the geometry optimization
    xtb_cmd = ["xtb", molecule_path_xyz, "--opt"]
    #result = subprocess.run(xtb_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    #attributes the outputs of the previous command
    #if result.returncode != 0:
     #   print("xTB failed:")
      #  print(result.stderr)
       # raise RuntimeError(f"xTB optimization failed:\n{result.stderr}")
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
    if not os.path.exists(output_xyz):
        raise FileNotFoundError(f"{output_xyz} not found after xTB optimization")

    
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



# Titles and site configuration
st.title('Create your reaction')
st.caption("Practical Proramming In Chemistry Project")
st.markdown("Draw two molecules and see if the product is stable enough")

# Exemple pour tester sur streamlit à enelever plus tard
#E_reagent1 = 150.0  
#E_reagent2 = 100.0  
#E_Prod = 200.0 

# Reset management
if "mol1" not in st.session_state:
    st.session_state.mol1 = ""
if "show_mol2" not in st.session_state:  
    st.session_state.show_mol2 = ""
if "mol2" not in st.session_state:
    st.session_state.mol2 = ""
if "product" not in st.session_state:
    st.session_state.product = ""



def draw_molecule(title, session_key):
    """Draw a molecule, show SMILES, 3D and energy."""
    st.subheader(f"Draw {title}") #titre
    mol_smiles = st_ketcher(st.session_state.get(session_key, ""), key=f"{session_key}_ketcher", height=400) #screen
    if mol_smiles:
        st.session_state[session_key] = mol_smiles
        with st.expander(f"SMILES {title}"):
            st.code(mol_smiles) # Smiles in a box
        mol_3D = generate_3D(mol_smiles)
        if mol_3D:
            with st.expander(f"3D Visualization {title}", expanded=False):
                visualize_3D(mol_3D) #3D molecule in a box
            try:
                elements, coords = smiles_to_3d(mol_smiles)
                xyz_filename=f"{session_key}.xyz"
                write_xyz_file(elements, coords, xyz_filename)
                with open(xyz_filename, "rb") as f:
                    st.download_button(
                        label=f"Download {title} .xyz file",
                        data=f,
                        file_name=xyz_filename,
                        mime="chemical/x-xyz"
                    )

                energy = GeomOptxyz_Energy(xyz_filename)
                st.session_state[f"energy_{session_key}"] = energy

            except Exception as e:
                st.error(f"Energy calculation failed: {str(e)}")
                energy = None

        with st.expander(f"Energy {title}"):
            if f"energy_{session_key}" in st.session_state:
                st.markdown(f"Energy: {st.session_state[f'energy_{session_key}']:.6f} Hartree")
            else:
                st.markdown("Energy not calculated yet")

    return mol_smiles




mol1 = draw_molecule("Molecule 1", "mol1") # draw molecule 1
if mol1:
    if st.button("Start to draw molecule 2"):
        st.session_state.show_mol2 = True # push the button to start molecule 2
    if st.session_state.show_mol2:
        mol2 = draw_molecule("Molecule 2", "mol2") # draw molecule 2


# Generated Product from Molecule 1 and Molecule 2
if st.session_state.mol1 and st.session_state.mol2:
    if st.button("Generate Product"):
        st.subheader("Generated Product from Molecule 1 and Molecule 2")
        product_smiles = st.session_state.mol1 + "." + st.session_state.mol2
        product = st_ketcher(product_smiles, key="product_ketcher", height=400)
        
        with st.expander("SMILES Product"):
            st.code(product_smiles)
            
        product_3D = generate_3D(product_smiles)
        if product_3D:
            with st.expander("Visualisation 3D Product", expanded=False):
                visualize_3D(product_3D)
                
            # Calculate product energy
            try:
                elements, coords = smiles_to_3d(product_smiles)
                temp_xyz = "temp_product.xyz"
                write_xyz_file(elements, coords, temp_xyz)
                E_Prod = GeomOptxyz_Energy(temp_xyz)
                os.remove(temp_xyz)
            except Exception as e:
                st.error(f"Product energy calculation failed: {str(e)}")
                E_Prod = None
        
        with st.expander("Energy Product"):
            if E_Prod is not None:
                st.markdown(f"Energy: {E_Prod:.6f} Hartree")
            else:
                st.markdown("Energy calculation failed")

        # Perform energy comparison if all energies are available
        if ('energy_mol1' in st.session_state and 
            'energy_mol2' in st.session_state and 
            E_Prod is not None):
            
            result = Energy_comparison(
                st.session_state['energy_mol1'],
                st.session_state['energy_mol2'],
                E_Prod
            )
            
            if result is not None:
                if result:
                    st.success("The reaction is thermodynamically stable at 0 K")
                else:
                    st.error("The reaction is not thermodynamically stable at 0 K")
            else:
                st.info("The reaction leads to a thermodynamic equilibrium at 0 K")
        else:
            st.warning("Could not perform energy comparison - missing energy values")

if st.button("Refresh Session State"):  
    st.session_state.mol1 = ""
    st.session_state.mol2 = ""
    st.session_state.show_mol2 = ""
    st.session_state.product = ""
    st.write("Session has been reset. Please refresh the page (F5).")
    
    































