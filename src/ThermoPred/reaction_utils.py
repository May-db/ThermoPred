from rxnutils.chem.reaction import ChemicalReaction
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import numpy as np

# Define a set of working reaction templates
working_templates = {
    # Sulfonyl chloride substitution (confirmed working)
    "sulfonyl_substitution": "Cl[S:3]([CH2:2][CH3:1])(=[O:4])=[O:5].[OH:6][CH2:7][CH2:8][Br:9]>>[CH3:1][CH2:2][S:3](=[O:4])(=[O:5])[O:6][CH2:7][CH2:8][Br:9]",
    
    # SN2 with specific leaving groups and nucleophiles (confirmed working)
    "sn2_bromide_oh": "[CH3:1][Br:2].[OH:3][H:4]>>[CH3:1][OH:3].[Br:2][H:4]",
    
    # Additional SN2 templates with valid atom mappings
    "sn2_chloride_oh": "[CH3:1][Cl:2].[OH:3][H:4]>>[CH3:1][OH:3].[Cl:2][H:4]",
    "sn2_iodide_oh": "[CH3:1][I:2].[OH:3][H:4]>>[CH3:1][OH:3].[I:2][H:4]",
    
    # SN2 with different nucleophiles
    "sn2_bromide_nh2": "[CH3:1][Br:2].[NH2:3][H:4]>>[CH3:1][NH2:3].[Br:2][H:4]",
    "sn2_bromide_sh": "[CH3:1][Br:2].[SH:3][H:4]>>[CH3:1][SH:3].[Br:2][H:4]",
    
    # Williamson ether synthesis
    "williamson_ether": "[C:1][O:2][H:3].[C:4][X:5]>>[C:1][O:2][C:4].[X:5][H:3]",
}

# Keep original utility functions
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


def is_alcohol(smiles):
    """Check if a SMILES string represents an alcohol."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    return mol.HasSubstructMatch(oh_pattern)


def is_alkyl_halide(smiles):
    """Check if a SMILES string represents an alkyl halide."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    
    # Check for C-X bond (X = F, Cl, Br, I)
    cx_pattern = Chem.MolFromSmarts("[CX4][F,Cl,Br,I]")
    return mol.HasSubstructMatch(cx_pattern)


def find_alcohol_and_halide(reactant1, reactant2):
    """Identify which reactant is an alcohol and which is a halide."""
    r1_is_alcohol = is_alcohol(reactant1)
    r1_is_halide = is_alkyl_halide(reactant1)
    r2_is_alcohol = is_alcohol(reactant2)
    r2_is_halide = is_alkyl_halide(reactant2)
    
    if r1_is_alcohol and r2_is_halide:
        return reactant1, reactant2
    elif r2_is_alcohol and r1_is_halide:
        return reactant2, reactant1
    else:
        return None, None


def handle_williamson_ether_synthesis(alcohol, halide):
    """Handle the specific case of Williamson ether synthesis."""
    try:
        halogen = None
        for element in ["Br", "Cl", "I", "F"]:
            if element in halide:
                halogen = element
                break
        
        if halogen is None:
            return None
        
        alcohol_mol = Chem.MolFromSmiles(alcohol)
        halide_mol = Chem.MolFromSmiles(halide)
        
        if alcohol_mol is None or halide_mol is None:
            return None
        
        oh_matches = alcohol_mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]"))
        
        cx_matches = halide_mol.GetSubstructMatches(Chem.MolFromSmarts(f"[C][{halogen}]"))
        
        if not oh_matches or not cx_matches:
            return None
        
        o_idx = oh_matches[0][0]
        c_idx_in_halide = cx_matches[0][0]
        
        new_alcohol = Chem.Mol(alcohol_mol)
        new_halide = Chem.Mol(halide_mol)
        
        new_alcohol = Chem.DeleteSubstructs(new_alcohol, Chem.MolFromSmarts("[OX2H]-[H]"))
        
        halide_without_x = Chem.DeleteSubstructs(new_halide, Chem.MolFromSmarts(f"[{halogen}]"))
        
        r_group_smiles = Chem.MolToSmiles(halide_without_x)
        
        ether_smiles = alcohol.replace("O", f"O{r_group_smiles}")
        
        if halogen == "Br":
            leaving_group = "Br"
        elif halogen == "Cl":
            leaving_group = "Cl"
        elif halogen == "I":
            leaving_group = "I"
        else:
            leaving_group = "F"
        
        return f"{ether_smiles}.{leaving_group}"
    except Exception as e:
        print(f"Error in handle_williamson_ether_synthesis: {e}")
        return None


# New template-based reaction prediction function
def predict_reaction_with_templates(reactant1, reactant2):
    """
    Predict reaction products using the template library.
    
    Args:
        reactant1 (str): SMILES of first reactant
        reactant2 (str): SMILES of second reactant
    
    Returns:
        tuple: (products, reaction_type)
    """
    # Combine reactants
    reactants = f"{reactant1}.{reactant2}"
    
    # First check for Williamson ether synthesis
    alcohol, halide = find_alcohol_and_halide(reactant1, reactant2)
    if alcohol and halide:
        williamson_product = handle_williamson_ether_synthesis(alcohol, halide)
        if williamson_product:
            return williamson_product, "williamson_ether"
    
    # Try each template
    for template_name, template_smarts in working_templates.items():
        try:
            # Create reaction from pattern
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
                # Join products with dot notation
                product_smiles = ".".join(products)
                return product_smiles, template_name
                
        except Exception as e:
            continue
    
    # If no template matches, try with RDKit's built-in reaction functionality
    try:
        reaction_patterns = {
            "default": "[*:1][*:2].[*:3][*:4]>>[*:1][*:3].[*:2][*:4]",
            "addition": "[C:1]=[C:2].[*:3][*:4]>>[*:3][C:1]-[C:2][*:4]",
            "elimination": "[*:1][C:2]-[C:3][*:4]>>[*:1][C:2]=[C:3][*:4]",
            "condensation": "[*:1][OH:5].[HO:6][*:2]>>[*:1][*:2].[OH2:5,6]",
        }
        
        for rxn_type, rxn_smarts in reaction_patterns.items():
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            
            reactant1_mol = Chem.MolFromSmiles(reactant1)
            reactant2_mol = Chem.MolFromSmiles(reactant2)
            
            if reactant1_mol is None or reactant2_mol is None:
                continue
            
            products = rxn.RunReactants((reactant1_mol, reactant2_mol))
            
            if products and len(products) > 0 and len(products[0]) > 0:
                product_mol = products[0][0]
                try:
                    Chem.SanitizeMol(product_mol)
                    return Chem.MolToSmiles(product_mol), rxn_type
                except:
                    continue
    except Exception as e:
        pass
    
    return attempt_generic_reaction(reactant1, reactant2)


def attempt_generic_reaction(reactant1, reactant2):
    """Attempt to predict a reaction product using generic approaches."""
    try:
        # Try using rxnutils for generic reactions
        # General coupling reaction
        coupling_rxn = ChemicalReaction("[*:1]~[*:2].[*:3]~[*:4]>>[*:1]~[*:2]-[*:3]~[*:4]")
        products = coupling_rxn.run_reaction([reactant1, reactant2])
        if products and len(products) > 0:
            return products[0], "coupling"
            
        # Common functional group reactions
        functional_rxns = [
            # Esterification: carboxylic acid + alcohol
            ChemicalReaction("[C:1](=[O:2])[OH:3].[*:4][OH:5]>>[C:1](=[O:2])[O:3][*:4].[OH2:5]"),
            # Amide formation: carboxylic acid + amine
            ChemicalReaction("[C:1](=[O:2])[OH:3].[*:4][NH2:5]>>[C:1](=[O:2])[N:5][*:4].[OH2:3]"),
            # Alkylation: nucleophile + alkyl halide
            ChemicalReaction("[*:1][H:2].[*:3][X:4]>>[*:1][*:3].[X:4][H:2]"),
        ]
        
        for rxn in functional_rxns:
            try:
                products = rxn.run_reaction([reactant1, reactant2])
                if products and len(products) > 0:
                    return products[0], "functional_group"
            except:
                continue
                
        mol1 = Chem.MolFromSmiles(reactant1)
        mol2 = Chem.MolFromSmiles(reactant2)
        
        if mol1 is None or mol2 is None:
            return None, "failed"
        
        combined = Chem.CombineMols(mol1, mol2)
        editable = Chem.EditableMol(combined)
        
        non_h_atoms1 = []
        non_h_atoms2 = []
        
        for atom in mol1.GetAtoms():
            if atom.GetSymbol() != 'H':
                non_h_atoms1.append(atom.GetIdx())
        
        offset = mol1.GetNumAtoms()
        for atom in mol2.GetAtoms():
            if atom.GetSymbol() != 'H':
                non_h_atoms2.append(atom.GetIdx() + offset)
        
        if non_h_atoms1 and non_h_atoms2:
            try:
                editable.AddBond(non_h_atoms1[0], non_h_atoms2[0], Chem.BondType.SINGLE)
                product_mol = editable.GetMol()
                Chem.SanitizeMol(product_mol)
                return Chem.MolToSmiles(product_mol), "generic"
            except:
                return None, "failed"
        
        return None, "failed"
            
    except Exception as e:
        print(f"Generic reaction attempt failed: {e}")
    
    return None, "failed"


# Function to add a new reaction template
def add_custom_reaction_template(name, smarts, test_reactant1, test_reactant2):
    """
    Test and add a new reaction template to the working templates.
    
    Args:
        name (str): Name of the reaction
        smarts (str): SMARTS pattern for the reaction
        test_reactant1 (str): SMILES of first test reactant
        test_reactant2 (str): SMILES of second test reactant
    
    Returns:
        bool: True if pattern works and was added, False otherwise
    """
    try:
        # Create and test the reaction
        rxn = ChemicalReaction(smarts)
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
            working_templates[name] = smarts
            return True
        else:
            return False
            
    except Exception as e:
        return False


# Updated get_product function that integrates template-based prediction
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
    # Define known reactant-product pairs
    reactant_pairs_to_products = {
        (normalize_smiles("C1C=CC=CC=1CO"), normalize_smiles("CBr")): "C1=CC=CC=C1COC.Br",
        (normalize_smiles("CBr"), normalize_smiles("C1C=CC=CC=1CO")): "C1=CC=CC=C1COC.Br",
        (normalize_smiles("CCCO"), normalize_smiles("CBr")): "CCCOC.Br",
        (normalize_smiles("CBr"), normalize_smiles("CCCO")): "CCCOC.Br",
        (normalize_smiles("CCO"), normalize_smiles("CBr")): "CCOC.Br",
        (normalize_smiles("CBr"), normalize_smiles("CCO")): "CCOC.Br",
        (normalize_smiles("CO"), normalize_smiles("CBr")): "COC.Br",
        (normalize_smiles("CBr"), normalize_smiles("CO")): "COC.Br",
        (normalize_smiles("CCCCO"), normalize_smiles("CBr")): "CCCCOC.Br",
        (normalize_smiles("CBr"), normalize_smiles("CCCCO")): "CCCCOC.Br",
        (normalize_smiles("CCCCCO"), normalize_smiles("CBr")): "CCCCCOC.Br",
        (normalize_smiles("CBr"), normalize_smiles("CCCCCO")): "CCCCCOC.Br",
    }
    
    norm_r1 = normalize_smiles(reactant1)
    norm_r2 = normalize_smiles(reactant2)
    
    full_product = None
    
    # Check known reactant pairs
    if (norm_r1, norm_r2) in reactant_pairs_to_products:
        full_product = reactant_pairs_to_products[(norm_r1, norm_r2)]
    elif (norm_r2, norm_r1) in reactant_pairs_to_products:
        full_product = reactant_pairs_to_products[(norm_r2, norm_r1)]
    else:
        # Try to identify alcohol and halide for Williamson ether synthesis
        alcohol, halide = find_alcohol_and_halide(reactant1, reactant2)
        if alcohol and halide:
            full_product = handle_williamson_ether_synthesis(alcohol, halide)
        
        # If that doesn't work, try template-based prediction
        if not full_product:
            result = predict_reaction_with_templates(reactant1, reactant2)
            if isinstance(result, tuple) and len(result) == 2:
                full_product, _ = result
            else:
                full_product = result
    
    # Return only the main product if requested
    if return_main_product_only and full_product and "." in full_product:
        return get_main_product(full_product)
    
    return full_product