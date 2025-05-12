import os
import numpy as np
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rxnutils.chem.reaction import ChemicalReaction
import subprocess


def generate_3D(smiles):
    """Generate 3D coordinates for a molecule from SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    molblock = Chem.MolToMolBlock(mol)
    return molblock


def smiles_to_3d(smiles, add_H=True, optimize=True):
    """Convert SMILES to 3D coordinates with optional H addition and optimization."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    if add_H:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if optimize:
        AllChem.MMFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coordinates = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return elements, coordinates


def write_xyz_file(elements, coordinates, filename):
    """Write molecule coordinates to XYZ file format."""
    os.makedirs(os.path.dirname(filename), exist_ok=True)  
    with open(filename, 'w') as f:
        f.write(f"{len(elements)}\n\n")
        for element, coord in zip(elements, coordinates):
            f.write(f"{element} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")


def GeomOptxyz_Energy(xyz_path):
    """Run xTB geometry optimization on an XYZ file and return the energy."""
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"{xyz_path} not found")

    xtb_cmd = ["xtb", xyz_path, "--opt"]
    result = subprocess.run(xtb_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        raise RuntimeError("xTB optimization failed")

    output_xyz = "xtbopt.xyz"
    output_log = "xtb.out"

    if not os.path.exists(output_xyz):
        raise FileNotFoundError("Optimization did not produce output xyz file")

    energy_hartree = None
    with open(output_log, "r") as f:
        for line in f:
            if "TOTAL ENERGY" in line:
                energy_hartree = float(line.strip().split()[-2])
                break

    return energy_hartree


def Energy_comparison(E1, E2, Ep):
    """Compare energies to determine if reaction is favorable."""
    delta_E = E1 + E2 - Ep
    if delta_E > 0:
        return "stable"
    elif delta_E == 0:
        return "equilibrium"
    else:
        return "unstable"


# New functions for reaction prediction without dataset dependency
def predict_reaction(reactant1_smiles, reactant2_smiles, reaction_type="default"):
    reaction_patterns = {
        "default": "[*:1][*:2].[*:3][*:4]>>[*:1][*:3].[*:2][*:4]",  # Simple substitution
        "addition": "[C:1]=[C:2].[*:3][*:4]>>[*:3][C:1]-[C:2][*:4]",  # Addition to alkene
        "elimination": "[*:1][C:2]-[C:3][*:4]>>[*:1][C:2]=[C:3][*:4]",  # Elimination
        "condensation": "[*:1][OH].[HO][*:2]>>[*:1][*:2]",  # Condensation with water elimination
        "substitution": "[C:1][X:2].[*:3]>>[C:1][*:3]",  # Nucleophilic substitution
    }
    
    rxn_smarts = reaction_patterns.get(reaction_type, reaction_patterns["default"])
    
    try:
        rxn = ChemicalReaction(rxn_smarts)
        products = rxn.run_reaction([reactant1_smiles, reactant2_smiles])
        if products and len(products) > 0:
            return products[0]
        return None
    except Exception as e:
        # Fallback to RDKit if rxnutils fails
        try:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            
            reactant1 = Chem.MolFromSmiles(reactant1_smiles)
            reactant2 = Chem.MolFromSmiles(reactant2_smiles)
            
            if reactant1 is None or reactant2 is None:
                return None
            
            products = rxn.RunReactants((reactant1, reactant2))
            
            # If we have products, return the SMILES of the first product
            if products and len(products) > 0 and len(products[0]) > 0:
                product_mol = products[0][0]
                try:
                    Chem.SanitizeMol(product_mol)
                    return Chem.MolToSmiles(product_mol)
                except:
                    return None
        except:
            pass
        return None


def predict_product_with_templates(reactant1, reactant2):
   
    reaction_types = ["default", "addition", "elimination", "condensation", "substitution"]
    
    for rxn_type in reaction_types:
        product = predict_reaction(reactant1, reactant2, rxn_type)
        if product:
            return product, rxn_type
    
    
    return attempt_generic_reaction(reactant1, reactant2)


def attempt_generic_reaction(reactant1, reactant2):
    
    try:
        # try using rxnutils for generic reactions
        # general coupling reaction
        coupling_rxn = ChemicalReaction("[*:1]~[*:2].[*:3]~[*:4]>>[*:1]~[*:2]-[*:3]~[*:4]")
        products = coupling_rxn.run_reaction([reactant1, reactant2])
        if products and len(products) > 0:
            return products[0], "coupling"
            
        # Common functional group reactions
        functional_rxns = [
            # Esterification: carboxylic acid + alcohol
            ChemicalReaction("[C:1](=[O:2])[OH].[*:3][OH]>>[C:1](=[O:2])[O][*:3]"),
            # Amide formation: carboxylic acid + amine
            ChemicalReaction("[C:1](=[O:2])[OH].[*:3][NH2]>>[C:1](=[O:2])[NH][*:3]"),
            # Alkylation: nucleophile + alkyl halide
            ChemicalReaction("[*:1][H].[*:2][Cl,Br,I]>>[*:1][*:2]"),
            # Reductive amination: aldehyde/ketone + amine
            ChemicalReaction("[C:1](=[O:2])[*:3].[NH2:4][*:5]>>[C:1]([NH:4][*:5])[*:3]"),
        ]
        
        for rxn in functional_rxns:
            try:
                products = rxn.run_reaction([reactant1, reactant2])
                if products and len(products) > 0:
                    return products[0], "functional_group"
            except:
                continue
                
        # If all else fails, fallback to RDKit approach
        # Create molecules from SMILES
        mol1 = Chem.MolFromSmiles(reactant1)
        mol2 = Chem.MolFromSmiles(reactant2)
        
        if mol1 is None or mol2 is None:
            return None, "failed"
        
        # Add atom mapping to distinguish atoms in the combined molecule
        for atom in mol1.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        
        offset = mol1.GetNumAtoms()
        for atom in mol2.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + offset + 1)
            
        # Try a simple bond formation between first atoms of each molecule
        combined = Chem.CombineMols(mol1, mol2)
        editable = Chem.EditableMol(combined)
        
        # Try to add a bond between first atoms of each molecule
        try:
            editable.AddBond(0, mol1.GetNumAtoms(), Chem.BondType.SINGLE)
            product_mol = editable.GetMol()
            # Remove atom mapping for clarity
            for atom in product_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            Chem.SanitizeMol(product_mol)
            return Chem.MolToSmiles(product_mol), "generic"
        except:
            return None, "failed"
            
    except Exception as e:
        print(f"Generic reaction attempt failed: {e}")
    
    return None, "failed"


# Legacy functions for compatibility
def load_reaction_data(csv_path):
    """function to maintain compatibility."""
    # This is just a placeholder to keep the interface consistent
    try:
        if os.path.exists(csv_path):
            return pd.read_csv(csv_path)
        else:
            return pd.DataFrame()  # Return empty DataFrame
    except:
        return pd.DataFrame()  # Return empty DataFrame


def get_product(df, reactant1, reactant2, rxn_column=None):
    """Get product from reactants using rxnutils and RDKit prediction methods."""
    # First check if we can use rxnutils to directly predict the reaction
    try:
        # Try a general reaction template with rxnutils
        general_rxn = ChemicalReaction("[*:1]~[*:2].[*:3]~[*:4]>>[*:1]~[*:3].[*:2]~[*:4]")
        products = general_rxn.run_reaction([reactant1, reactant2])
        if products and len(products) > 0:
            return products[0]
    except:
        pass
    
    # If that fails rdkit method
    product, rxn_type = predict_product_with_templates(reactant1, reactant2)
    return product  # Return only the product SMILES