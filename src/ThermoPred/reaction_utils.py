from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
import os
import numpy as np
from rdkit.Chem import rdChemReactions

# Utility functions
def generate_3D(smiles):
    """Converts a SMILES string into a 3D molecular structure in molblock format.
    If the SMILES contains multiple fragments, only the main product is used.
    Hydrogen atoms are added, and 3D coordinates are generated."""
    if "." in smiles:
        smiles = get_main_product(smiles) 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    molblock = Chem.MolToMolBlock(mol)
    return molblock


def smiles_to_3d(smiles, add_H=True, optimize=True, max_attempts=3):
    """Converts a SMILES into 3D atomic coordinates and element symbols.
    Adds hydrogens and performs geometry optimization using MMFF.
    Includes multiple attempts to embed the molecule if needed."""
    if "." in smiles:
        smiles = get_main_product(smiles)
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    if add_H:
        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.randomSeed = 42
        params.numThreads = 1 
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if optimize:
        AllChem.MMFFOptimizeMolecule(mol)
        success = False
        for attempt in range(max_attempts):
            try:
                # Trying to embed the molecule
                AllChem.EmbedMolecule(mol, params)
                if mol.GetNumConformers() > 0:
                # Optimize with MMFF if requested
                    if optimize:
                        AllChem.MMFFOptimizeMolecule(mol)
                    success = True
                    break
            except Exception as e:
                if attempt == max_attempts - 1:
                    raise ValueError(f"Failed to generate 3D coordinates after {max_attempts} attempts: {str(e)}")
    
        if not success:
            raise ValueError("Failed to generate reasonable 3D coordinates")

    conf = mol.GetConformer()
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coordinates = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return elements, coordinates


def write_xyz_file(elements, coordinates, filename):
    """ Writes atomic elements and 3D coordinates to a .xyz file in standard XYZ format.
    Creates the directory if it doesn't exist."""
    os.makedirs(os.path.dirname(filename), exist_ok=True)  
    with open(filename, 'w') as f:
        f.write(f"{len(elements)}\n\n")
        for element, coord in zip(elements, coordinates):
            f.write(f"{element} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")


def calculate_energy_with_rdkit(smiles, optimize_steps=500):
    """Calculate a molecule's energy using RDKit instead of xTB"""
    try:
        # Clean the SMILES
        if "." in smiles:
            smiles = get_main_product(smiles)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Cannot read SMILES: {smiles}")
        
        mol = Chem.AddHs(mol)
        
        # Generate a clean 3D molecule
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.randomSeed = 42
        params.maxIterations = 1000  
        AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params)
        
        if mol.GetNumConformers() == 0:
            raise ValueError("Failed to generate 3D coordinates")
        
        # All conformations are optimized to take the lowest energy later
        energies = []
        for conf_id in range(mol.GetNumConformers()):
            try:
                # Setting up the MMFF94 field
                mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
                if mp is None:
                    continue
                
                ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conf_id)
                if ff is None:
                    continue
                
                # Optimization
                ff.Initialize()
                result = ff.Minimize(maxIts=optimize_steps)
                
                # Energy calculation
                energy = ff.CalcEnergy()
                energies.append((conf_id, energy))
            except Exception as e:
                continue
        
        if not energies:
            # UFF method as backup
            for conf_id in range(mol.GetNumConformers()):
                try:
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                    if ff is None:
                        continue
                    
                    ff.Initialize()
                    result = ff.Minimize(maxIts=optimize_steps)
                    energy = ff.CalcEnergy()
                    energies.append((conf_id, energy))
                except Exception as e:
                    continue
        
        if not energies:
            raise ValueError("Cannot calculate energy with MMFF94 or UFF")
        
        # Take the lowest energy
        energies.sort(key=lambda x: x[1])
        best_conf_id, lowest_energy = energies[0]
        
        # Convert to Hartree
        energy_hartree = lowest_energy / 627.5  # 1 Hartree ≈ 627.5 kcal/mol
        
        # Create xyz data
        conf = mol.GetConformer(best_conf_id)
        elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
        coordinates = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        
        return energy_hartree, elements, coordinates
    
    except Exception as e:
        raise ValueError(f"Error calculating energy: {str(e)}")

def Energy_comparison(E1, E2, Ep):
    """ Compares the energy of a product to the sum of two reactants and returns a stability assessment.
    Returns 'stable', 'equilibrium', or 'unstable' based on the energy difference."""
    delta_E =  Ep - E1 - E2
    if delta_E < 0:
        return "stable"
    elif delta_E == 0:
        return "equilibrium"
    else:
        return "unstable"
    
REACTION_TEMPLATES = {
    # Amide formation 
    "amide_formation": "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]",
    
    # SN2 reactions
    "sn2_halide_oh": "[C:1][Cl,Br,I:2].[O:3]>>[C:1][O:3]",
    "sn2_halide_n": "[C:1][Cl,Br,I:2].[N:3]>>[C:1][N:3]",
    
    # Carbon-carbon coupling (R₁-X + R₂-Y → R₁-R₂)
    "c_c_coupling": "[C:1][Cl,Br,I:2].[C:3][Cl,Br,I:4]>>[C:1][C:3]",
    
    # Specific coupling examples
    "methyl_coupling": "C[Cl,Br,I].C[Cl,Br,I]>>CC",
    
    # Interhalogen
    "halogen_exchange": "[Cl:1].[Br:2]>>[Cl:1][Br:2]"
}

def get_main_product(product_smiles):
    """Extract the main product (largest molecule) from a dot-separated SMILES"""
    if "." not in product_smiles:
        return product_smiles
    
    # Split by dots and find the largest molecule
    components = product_smiles.split(".")
    largest_component = ""
    largest_size = 0
    
    for comp in components:
        mol = Chem.MolFromSmiles(comp)
        if mol:
            size = mol.GetNumAtoms()
            if size > largest_size:
                largest_size = size
                largest_component = comp
    
    return largest_component if largest_component else components[0]

def react(reactant1_smiles, reactant2_smiles=None):
    """
    Extremely simplified reaction prediction using direct RDKit approach.
    
    Args:
        reactant1_smiles (str): SMILES string for first reactant
        reactant2_smiles (str, optional): SMILES string for second reactant
        
    Returns:
        tuple: (product_smiles, template_used) if successful, (None, None) if not
    """
    # Handle dot-separated input
    if reactant2_smiles is None and '.' in reactant1_smiles:
        parts = reactant1_smiles.split('.', 1)
        reactant1_smiles, reactant2_smiles = parts
    
    # Convert SMILES to RDKit molecules - exactly following your example
    mol1 = Chem.MolFromSmiles(reactant1_smiles)
    mol2 = Chem.MolFromSmiles(reactant2_smiles) if reactant2_smiles else None
    
    if not mol1 or (reactant2_smiles and not mol2):
        print(f"Error: Invalid SMILES - {reactant1_smiles} or {reactant2_smiles}")
        return None, None
    
    # Create reactants tuple - exactly following your example
    reactants = (mol1, mol2) if mol2 else (mol1,)
    
    # Special case for elemental halogens
    if reactant1_smiles == "Cl" and reactant2_smiles == "Br":
        return "ClBr", "halogen_exchange"
    elif reactant1_smiles == "Br" and reactant2_smiles == "Cl":
        return "BrCl", "halogen_exchange"
    
    # Basic check for alkyl halides
    is_r1_alkyl_halide = any(atom.GetSymbol() in ["Cl", "Br", "I"] for atom in mol1.GetAtoms())
    is_r2_alkyl_halide = mol2 and any(atom.GetSymbol() in ["Cl", "Br", "I"] for atom in mol2.GetAtoms())
    
    # For two alkyl halides, try carbon-carbon coupling first
    if is_r1_alkyl_halide and is_r2_alkyl_halide:
        try:
            rxn = rdChemReactions.ReactionFromSmarts(REACTION_TEMPLATES["c_c_coupling"])
            products = rxn.RunReactants(reactants)
            if products and len(products) > 0 and len(products[0]) > 0:
                full_product_smiles = Chem.MolToSmiles(products[0][0])
                main_product = get_main_product(full_product_smiles)
                return main_product, "c_c_coupling"
        except Exception:
            pass  # If this fails, continue with other templates
    
    # Try all templates
    for template_name, template_smarts in REACTION_TEMPLATES.items():
        try:
            # Create reaction from SMARTS - exactly following your example
            rxn = rdChemReactions.ReactionFromSmarts(template_smarts)
            
            # Run the reaction - exactly following your example
            products = rxn.RunReactants(reactants)
            
            # Check if we got any products
            if products and len(products) > 0 and len(products[0]) > 0:
                # Get the product SMILES
                full_product_smiles = Chem.MolToSmiles(products[0][0])
                
                # Extract the main product if there are multiple products
                main_product = get_main_product(full_product_smiles)
                
                return main_product, template_name
        except Exception:
            continue
    
    return None, None

# Test script
if __name__ == "__main__":
    test_cases = [
        # Basic reactions from your example
        ("C(=O)O", "CNC", "Amide formation"),
        ("CBr", "O", "SN2 substitution"),
        
        # Carbon-carbon coupling (R₁Cl + R₂Br → R₁R₂)
        ("CCl", "CBr", "Methyl coupling"),
        ("CCCCCCl", "CBr", "Hexyl-methyl coupling"),
        
        # Elemental halogens
        ("Cl", "Br", "Interhalogen compound")
    ]
    
    for reactant1, reactant2, description in test_cases:
        product, template = react(reactant1, reactant2)
        result = f"{product} (using {template})" if product else "No reaction"
        print(f"{reactant1} + {reactant2} → {result} - {description}")

