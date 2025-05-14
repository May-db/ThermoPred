from rxnutils.chem.reaction import ChemicalReaction
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import numpy as np

# Define a set of reaction patterns that work
working_templates = {
    # Sulfonyl chloride substitution (confirmed working)
    "sulfonyl_substitution": "Cl[S:3]([CH2:2][CH3:1])(=[O:4])=[O:5].[OH:6][CH2:7][CH2:8][Br:9]>>[CH3:1][CH2:2][S:3](=[O:4])(=[O:5])[O:6][CH2:7][CH2:8][Br:9]",
    
    # SN2 with specific leaving groups and nucleophiles (confirmed working)
    "sn2_bromide_oh": "[CH3:1][Br:2].[OH:3][H:4]>>[CH3:1][OH:3].[Br:2][H:4]",
    "sn2_bromide_oh_2": "[C:1][Br:2].[O:3][H:4]>>[C:1][O:3].[Br:2][H:4]",  # More general version
    
    # Additional SN2 templates with valid atom mappings
    "sn2_chloride_oh": "[CH3:1][Cl:2].[OH:3][H:4]>>[CH3:1][OH:3].[Cl:2][H:4]",
    "sn2_chloride_oh_2": "[C:1][Cl:2].[O:3][H:4]>>[C:1][O:3].[Cl:2][H:4]",  # More general version
    "sn2_iodide_oh": "[CH3:1][I:2].[OH:3][H:4]>>[CH3:1][OH:3].[I:2][H:4]",
    "sn2_iodide_oh_2": "[C:1][I:2].[O:3][H:4]>>[C:1][O:3].[I:2][H:4]",  # More general version
    
    # SN2 with different nucleophiles
    "sn2_bromide_nh2": "[CH3:1][Br:2].[NH2:3][H:4]>>[CH3:1][NH2:3].[Br:2][H:4]",
    "sn2_bromide_sh": "[CH3:1][Br:2].[SH:3][H:4]>>[CH3:1][SH:3].[Br:2][H:4]",
    
    # Alkyl halide coupling
    "alkyl_coupling": "[C:1][Cl:2].[C:3][Br:4]>>[C:1][C:3].[Cl:2][Br:4]"
}

# Utility functions
def generate_3D(smiles):
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
            smiles = smiles.split(".")[0]  # Take only the first molecule
        
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
    delta_E = E1 + E2 - Ep
    if delta_E < 0:
        return "stable"
    elif delta_E == 0:
        return "equilibrium"
    else:
        return "unstable"


def normalize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol)
    except:
        return smiles


def get_main_product(product_smiles):
    if product_smiles is None:
        return None
        
    if "." not in product_smiles:
        return product_smiles
    
    components = product_smiles.split(".")
    
    if len(components) == 1:
        return components[0]
    largest_component = components[0]
    largest_size = 0
    
    for component in components:
        mol = Chem.MolFromSmiles(component)
        if mol:
            size = mol.GetNumAtoms()
            if size > largest_size:
                largest_size = size
                largest_component = component
    
    return largest_component


def is_halide(smiles):
    """Check if a molecule contains a halide (F, Cl, Br, I)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I']:
            return True
    
    return False


def is_alkyl_halide(smiles):
    """Check if a SMILES string represents an alkyl halide."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    
    # Check for C-X bond (X = F, Cl, Br, I)
    cx_pattern = Chem.MolFromSmarts("[CX4][F,Cl,Br,I]")
    return mol.HasSubstructMatch(cx_pattern)


def is_alcohol(smiles):
    """Check if a SMILES string represents an alcohol."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    return mol.HasSubstructMatch(oh_pattern)


def predict_sn2_with_hydroxide(alkyl_halide, hydroxide):
    """Specialized function to predict SN2 reactions with hydroxide."""
    # Find the halogen
    halogen = None
    for h in ['F', 'Cl', 'Br', 'I']:
        if h in alkyl_halide:
            halogen = h
            break
    
    if halogen:
        # Replace the halogen with OH
        alkyl_group = alkyl_halide.replace(halogen, "")
        alcohol_product = alkyl_group + "O"
        
        # Check if the product is valid
        mol = Chem.MolFromSmiles(alcohol_product)
        if mol:
            return f"{alcohol_product}.{halogen}"
    
    return None


def predict_reaction_products(reactant1, reactant2):
    """
    Predict the products of a reaction between two reactants using working templates.
    
    Args:
        reactant1 (str): SMILES of first reactant
        reactant2 (str): SMILES of second reactant
    
    Returns:
        dict: Dictionary of {reaction_type: predicted_products}
    """
    # Check if one is an alkyl halide and the other is hydroxide
    r1_is_alkyl_halide = is_alkyl_halide(reactant1)
    r2_is_alkyl_halide = is_alkyl_halide(reactant2)
    
    # Simple check for hydroxide
    r1_is_hydroxide = reactant1 == "O" or reactant1 == "[OH-]" or reactant1 == "OH"
    r2_is_hydroxide = reactant2 == "O" or reactant2 == "[OH-]" or reactant2 == "OH"
    
    # SN2 reaction with hydroxide
    if (r1_is_alkyl_halide and r2_is_hydroxide):
        result = predict_sn2_with_hydroxide(reactant1, reactant2)
        if result:
            return {"sn2_hydroxide": [result]}
    elif (r2_is_alkyl_halide and r1_is_hydroxide):
        result = predict_sn2_with_hydroxide(reactant2, reactant1)
        if result:
            return {"sn2_hydroxide": [result]}
    
    # Check for alkyl halide coupling
    if r1_is_alkyl_halide and r2_is_alkyl_halide:
        # Get the carbon count
        carbon_count_r1 = reactant1.count('C')
        carbon_count_r2 = reactant2.count('C')
        
        # Find halogens
        halogen_r1 = next((h for h in ['F', 'Cl', 'Br', 'I'] if h in reactant1), None)
        halogen_r2 = next((h for h in ['F', 'Cl', 'Br', 'I'] if h in reactant2), None)
        
        if halogen_r1 and halogen_r2:
            # Create carbon chain product
            product = 'C' * (carbon_count_r1 + carbon_count_r2)
            
            # Create halogen product
            halogen_product = f"{halogen_r1}{halogen_r2}"
            
            results = {"alkyl_coupling": [f"{product}.{halogen_product}"]}
            return results
    
    # Combine reactants
    reactants = f"{reactant1}.{reactant2}"
    
    results = {}
    
    # Try each template
    for template_name, template_smarts in working_templates.items():
        try:
            rxn = ChemicalReaction(template_smarts)
            rxn.generate_reaction_template(radius=1)
            
            # Apply the template to predict products
            product_list = rxn.canonical_template.apply(reactants)
            
            # Flatten the product list
            products = []
            if product_list:
                for product_set in product_list:
                    products.extend(product_set)
            
            if products:
                results[template_name] = products
                
        except Exception as e:
            continue
    
    # If no results yet, try direct RDKit reaction approach
    if not results:
        try:
            rxn = AllChem.ReactionFromSmarts("[C:1][F,Cl,Br,I:2].[O:3][H:4]>>[C:1][O:3].[*:2][H:4]")
            r1_mol = Chem.MolFromSmiles(reactant1)
            r2_mol = Chem.MolFromSmiles(reactant2)
            
            if r1_mol and r2_mol:
                # Try both directions
                products = rxn.RunReactants((r1_mol, r2_mol))
                
                if not products or len(products) == 0 or len(products[0]) == 0:
                    products = rxn.RunReactants((r2_mol, r1_mol))
                
                if products and len(products) > 0 and len(products[0]) > 0:
                    product_smiles = []
                    for product_set in products:
                        valid_products = []
                        for p in product_set:
                            try:
                                Chem.SanitizeMol(p)
                                valid_products.append(Chem.MolToSmiles(p))
                            except:
                                pass
                        if valid_products:
                            product_smiles.append(".".join(valid_products))
                    
                    if product_smiles:
                        results["rdkit_sn2"] = product_smiles
        except:
            pass
    
    return results


def add_custom_sn2_pattern(r_group, leaving_group, nucleophile, test_reactant1, test_reactant2):
    """
    Add a custom SN2 pattern based on specific R-group, leaving group, and nucleophile.
    
    Args:
        r_group (str): R-group atom mapping (e.g., "[CH3:1]", "[CH2:1][CH3:5]")
        leaving_group (str): Leaving group atom mapping (e.g., "[Br:2]")
        nucleophile (str): Nucleophile atom mapping (e.g., "[OH:3]")
        test_reactant1 (str): Test reactant SMILES
        test_reactant2 (str): Test reactant SMILES
    
    Returns:
        bool: True if successful, False otherwise
    """
    # Create pattern name
    r_group_str = r_group.replace("[", "").replace("]", "").replace(":", "").split(")")[0]
    lg_str = leaving_group.replace("[", "").replace("]", "").replace(":", "").split(")")[0]
    nu_str = nucleophile.replace("[", "").replace("]", "").replace(":", "").split(")")[0]
    pattern_name = f"sn2_{r_group_str}_{lg_str}_{nu_str}"
    
    # Create SMARTS pattern
    pattern_smarts = f"{r_group}{leaving_group}.[{nucleophile}][H:4]>>{r_group}{nucleophile}.[{leaving_group}][H:4]"
    
    try:
        # Create and test the reaction
        rxn = ChemicalReaction(pattern_smarts)
        rxn.generate_reaction_template(radius=1)
        
        # Test with provided reactants
        reactants = f"{test_reactant1}.{test_reactant2}"
        
        # Apply the template
        product_list = rxn.canonical_template.apply(reactants)
        
        # Check results
        products = []
        if product_list:
            for product_set in product_list:
                products.extend(product_set)
        
        if products:
            working_templates[pattern_name] = pattern_smarts
            return True
        else:
            return False
            
    except Exception as e:
        return False


def get_product(reactant1, reactant2, rxn_column=None, return_main_product_only=False):
    """
    Predict the product of a reaction between two reactants.
    
    Args:
        reactant1 (str): SMILES of first reactant
        reactant2 (str): SMILES of second reactant
        rxn_column (str, optional): Not used, kept for compatibility
        return_main_product_only (bool): Whether to return only the main product
    
    Returns:
        str: SMILES of the predicted product
    """
    # Check if one is an alkyl halide and the other is hydroxide
    r1_is_alkyl_halide = is_alkyl_halide(reactant1)
    r2_is_alkyl_halide = is_alkyl_halide(reactant2)
    
    # Simple check for hydroxide
    r1_is_hydroxide = reactant1 == "O" or reactant1 == "[OH-]" or reactant1 == "OH"
    r2_is_hydroxide = reactant2 == "O" or reactant2 == "[OH-]" or reactant2 == "OH"
    
    # SN2 reaction with hydroxide
    if (r1_is_alkyl_halide and r2_is_hydroxide):
        result = predict_sn2_with_hydroxide(reactant1, reactant2)
        if result:
            if return_main_product_only:
                return result.split(".")[0]
            return result
    elif (r2_is_alkyl_halide and r1_is_hydroxide):
        result = predict_sn2_with_hydroxide(reactant2, reactant1)
        if result:
            if return_main_product_only:
                return result.split(".")[0]
            return result
    
    # Check for alkyl halide coupling first
    if r1_is_alkyl_halide and r2_is_alkyl_halide:
        # Get the carbon count
        carbon_count_r1 = reactant1.count('C')
        carbon_count_r2 = reactant2.count('C')
        
        # Find halogens
        halogen_r1 = next((h for h in ['F', 'Cl', 'Br', 'I'] if h in reactant1), None)
        halogen_r2 = next((h for h in ['F', 'Cl', 'Br', 'I'] if h in reactant2), None)
        
        if halogen_r1 and halogen_r2:
            # Create carbon chain product
            product = 'C' * (carbon_count_r1 + carbon_count_r2)
            
            # Create halogen product
            halogen_product = f"{halogen_r1}{halogen_r2}"
            
            full_product = f"{product}.{halogen_product}"
            
            if return_main_product_only:
                return product
            else:
                return full_product
    
    # Get products using template-based approach
    results = predict_reaction_products(reactant1, reactant2)
    
    # If no templates matched, return None
    if not results:
        return None
    
    # Take the first matching template's products
    first_template = next(iter(results))
    products = results[first_template]
    
    # Handle different return types
    if isinstance(products, list):
        # Join products with dot notation
        full_product = ".".join(products) if len(products) > 0 else None
    else:
        full_product = products
    
    # Return only the main product if requested
    if return_main_product_only and full_product and "." in full_product:
        return get_main_product(full_product)
    
    return full_product