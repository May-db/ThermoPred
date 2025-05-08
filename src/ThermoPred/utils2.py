from morfeus import read_xyz
import subprocess
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import streamlit as st
from stmol import showmol
import py3Dmol
from streamlit_ketcher import st_ketcher

# Constants
CALCULATION_TIMEOUT = 120  # seconds

def initialize_session_state():
    """Initialize all required session state variables"""
    defaults = {
        "mol1": "",
        "mol2": "",
        "product": "",
        "show_mol2": False,
        "energy_mol1": None,
        "energy_mol2": None,
        "energy_product": None
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

def generate_3d_structure(smiles):
    """Generate 3D coordinates from SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
            
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        AllChem.EmbedMolecule(mol, params)
        
        return Chem.MolToMolBlock(mol)
    except Exception as e:
        st.error(f"3D structure generation failed: {str(e)}")
        return None

def visualize_molecule(mol_block):
    """Visualize molecule in 3D using py3Dmol"""
    try:
        view = py3Dmol.view(width=400, height=400)
        view.addModel(mol_block, 'mol')
        view.setStyle({
            'sphere': {'colorscheme': 'cyanCarbon', 'scale': 0.25},
            'stick': {'colorscheme': 'cyanCarbon'}
        })
        view.zoomTo()
        view.spin()
        view.setBackgroundColor('white')
        showmol(view, height=400, width=400)
    except Exception as e:
        st.error(f"Visualization failed: {str(e)}")

def create_xyz_file(elements, coordinates):
    """Create temporary XYZ file and return its path"""
    try:
        fd, path = tempfile.mkstemp(suffix='.xyz')
        with os.fdopen(fd, 'w') as f:
            f.write(f"{len(elements)}\n\n")
            for element, coord in zip(elements, coordinates):
                f.write(f"{element} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
        return path
    except Exception as e:
        raise RuntimeError(f"XYZ file creation failed: {str(e)}")

def run_xtb_calculation(xyz_path):
    """Run xTB calculation and return energy"""
    try:
        # Verify input file
        if not os.path.exists(xyz_path):
            raise FileNotFoundError(f"Input file not found: {xyz_path}")

        # Run xTB
        cmd = ["xtb", xyz_path, "--opt", "--chrg", "0", "--uhf", "0"]
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=CALCULATION_TIMEOUT
        )

        # Check results
        if result.returncode != 0:
            error_msg = [
                f"xTB failed (return code {result.returncode})",
                "=== Error Output ===",
                result.stderr,
                "=== Standard Output ===",
                result.stdout
            ]
            raise RuntimeError("\n".join(error_msg))

        # Parse energy from output
        output_log = "xtb.out"
        if not os.path.exists(output_log):
            raise FileNotFoundError("xTB output log not found")

        with open(output_log) as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    return float(line.split()[3])

        raise ValueError("Energy not found in xTB output")
    finally:
        # Clean up temporary files
        for f in ["xtbopt.xyz", "xtb.out", "xtbopt.log", "xtbrestart"]:
            if os.path.exists(f):
                os.remove(f)

def calculate_molecule_energy(smiles):
    """Calculate energy for a molecule given its SMILES"""
    try:
        # Generate 3D structure
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Get elements and coordinates
        conf = mol.GetConformer()
        elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
        coordinates = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        
        # Create temporary XYZ file
        xyz_path = create_xyz_file(elements, coordinates)
        
        # Run xTB calculation
        energy = run_xtb_calculation(xyz_path)
        
        # Clean up
        os.remove(xyz_path)
        
        return energy
    except Exception as e:
        raise RuntimeError(f"Energy calculation failed: {str(e)}")

def energy_comparison(E_reagent1, E_reagent2, E_product):
    """Compare energies and determine thermodynamic stability"""
    try:
        E_start = E_reagent1 + E_reagent2
        dE = E_start - E_product
        
        if dE > 0:
            return True, "The reaction is thermodynamically stable at 0 K"
        elif dE == 0:
            return None, "The reaction leads to a thermodynamic equilibrium at 0 K"
        else:
            return False, "The reaction is not thermodynamically stable at 0 K"
    except Exception as e:
        raise RuntimeError(f"Energy comparison failed: {str(e)}")

def draw_molecule_interface(title, session_key):
    """Create interface for drawing and analyzing a molecule"""
    st.subheader(f"Draw {title}")
    
    # Molecule drawing interface
    smiles = st_ketcher(
        st.session_state.get(session_key, ""),
        key=f"{session_key}_ketcher",
        height=400
    )
    
    if smiles:
        st.session_state[session_key] = smiles
        
        # Display SMILES
        with st.expander(f"SMILES {title}"):
            st.code(smiles)
        
        # Generate and display 3D structure
        mol_block = generate_3d_structure(smiles)
        if mol_block:
            with st.expander(f"3D Visualization {title}", expanded=False):
                visualize_molecule(mol_block)
            
            # Calculate and display energy
            try:
                energy = calculate_molecule_energy(smiles)
                st.session_state[f"energy_{session_key}"] = energy
                
                with st.expander(f"Energy {title}"):
                    st.markdown(f"Energy: {energy:.6f} Hartree")
            except Exception as e:
                st.error(f"Energy calculation failed: {str(e)}")
        
    return smiles

def main():
    """Main application function"""
    st.title('Reaction Thermodynamics Analyzer')
    st.caption("Practical Programming in Chemistry Project")
    st.markdown("Draw two molecules and analyze their reaction thermodynamics")
    
    # Initialize session state
    initialize_session_state()
    
    # Molecule 1 interface
    mol1 = draw_molecule_interface("Molecule 1", "mol1")
    
    # Molecule 2 interface (only shown after Molecule 1 is drawn)
    if mol1 and st.button("Start drawing Molecule 2"):
        st.session_state.show_mol2 = True
    
    if st.session_state.show_mol2:
        mol2 = draw_molecule_interface("Molecule 2", "mol2")
    
    # Reaction analysis (only shown when both molecules are drawn)
    if st.session_state.mol1 and st.session_state.mol2:
        if st.button("Analyze Reaction"):
            st.subheader("Reaction Analysis")
            
            # Generate product SMILES (simple concatenation for demo)
            product_smiles = f"{st.session_state.mol1}.{st.session_state.mol2}"
            st.session_state.product = product_smiles
            
            # Display product information
            with st.expander("Product SMILES"):
                st.code(product_smiles)
            
            # Generate and display product 3D structure
            product_3d = generate_3d_structure(product_smiles)
            if product_3d:
                with st.expander("Product 3D Visualization", expanded=False):
                    visualize_molecule(product_3d)
                
                # Calculate product energy
                try:
                    product_energy = calculate_molecule_energy(product_smiles)
                    st.session_state.energy_product = product_energy
                    
                    with st.expander("Product Energy"):
                        st.markdown(f"Energy: {product_energy:.6f} Hartree")
                except Exception as e:
                    st.error(f"Product energy calculation failed: {str(e)}")
            
            # Perform energy comparison if all energies are available
            if (st.session_state.get("energy_mol1") is not None and
                st.session_state.get("energy_mol2") is not None and
                st.session_state.get("energy_product") is not None):
                
                try:
                    is_stable, message = energy_comparison(
                        st.session_state.energy_mol1,
                        st.session_state.energy_mol2,
                        st.session_state.energy_product
                    )
                    
                    if is_stable is True:
                        st.success(message)
                    elif is_stable is False:
                        st.error(message)
                    else:
                        st.info(message)
                except Exception as e:
                    st.error(f"Energy comparison failed: {str(e)}")
            else:
                st.warning("Could not perform energy comparison - missing energy values")
    
    # Reset button
    if st.button("Reset All"):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.experimental_rerun()

if __name__ == "__main__":
    main()