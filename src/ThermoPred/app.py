# app.py
import os
import streamlit as st
from stmol import showmol
from streamlit_ketcher import st_ketcher
import py3Dmol
from pathlib import Path
from rdkit import Chem

# Import functions from the updated reaction_utils
from reaction_utils import (
    generate_3D, 
    smiles_to_3d, 
    write_xyz_file, 
    calculate_energy_with_rdkit, 
    Energy_comparison, 
    predict_reaction_products,
    working_templates,
    get_main_product,
    is_halide,
    is_alkyl_halide,
    is_alcohol
)

# Set up the Streamlit page
st.set_page_config(page_title="Reaction Thermo Tool", layout="centered")
st.title("🧪 Reaction Thermodynamics Predictor")
st.caption("Draw two molecules to predict if the reaction is thermodynamically favorable at 0 K.")

# Create directory for xyz files if it doesn't exist
os.makedirs("xyz_files", exist_ok=True)

def visualize_3D(molblock):
    """Create a 3D visualization of a molecule."""
    view = py3Dmol.view(width=400, height=400)
    view.addModel(molblock, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.25}})
    view.zoomTo()
    view.spin()
    view.setBackgroundColor('white')
    showmol(view, height=400, width=400)

def draw_and_process(title, session_key):
    """Draw a molecule and calculate its energy."""
    st.subheader(f"{title}")
    smiles = st_ketcher(st.session_state.get(session_key, ""), key=f"{session_key}_ketcher", height=400)
    
    if smiles:
        st.session_state[session_key] = smiles
        with st.expander(f"SMILES"):
            st.code(smiles)
        
        molblock = generate_3D(smiles)
        if molblock:
            with st.expander("3D Visualization"):
                visualize_3D(molblock)
            
            try:
                energy, elements, coords = calculate_energy_with_rdkit(smiles)
                
                # register the optimized geometry
                os.makedirs("xyz_files", exist_ok=True)
                xyz_path = os.path.abspath(f"xyz_files/{session_key}_optimized.xyz")
                write_xyz_file(elements, coords, xyz_path)
                
                st.session_state[f"{session_key}_energy"] = energy
                with st.expander("Energy"):
                    st.success(f"Energy (RDKit): {energy:.6f} Hartree")
                
            except Exception as e:
                st.error(f"Energy calculation failed: {str(e)}")
                import traceback
                st.code(traceback.format_exc())
    
    return smiles


if "mol1" not in st.session_state:
    st.session_state.mol1 = ""
if "mol2" not in st.session_state:
    st.session_state.mol2 = ""



mol1 = draw_and_process("Molecule 1", "mol1")
mol2 = draw_and_process("Molecule 2", "mol2")

# Display molecule types if identified
if mol1 and mol2:
    st.subheader("Reaction")
    types_expander = st.expander("Molecule Types", expanded=False)
    with types_expander:
        col1, col2 = st.columns(2)
        with col1:
            st.write("**Molecule 1:**")
            mol_types = []
            if is_alkyl_halide(mol1):
                mol_types.append("Alkyl halide (C-X)")
            elif is_halide(mol1):
                mol_types.append("Contains halide (F, Cl, Br, I)")
            if is_alcohol(mol1):
                mol_types.append("Contains alcohol group (-OH)")
            
            if mol_types:
                for t in mol_types:
                    st.write(f"- {t}")
            else:
                st.write("- No specific functional groups identified")
                
        with col2:
            st.write("**Molecule 2:**")
            mol_types = []
            if is_alkyl_halide(mol2):
                mol_types.append("Alkyl halide (C-X)")
            elif is_halide(mol2):
                mol_types.append("Contains halide (F, Cl, Br, I)")
            if is_alcohol(mol2):
                mol_types.append("Contains alcohol group (-OH)")
            
            if mol_types:
                for t in mol_types:
                    st.write(f"- {t}")
            else:
                st.write("- No specific functional groups identified")
                
        # Suggest possible reactions
        st.write("**Possible Reactions:**")
        if is_alkyl_halide(mol1) and is_alkyl_halide(mol2):
            st.write("- Alkyl-alkyl coupling (Wurtz reaction)")
        elif (is_alkyl_halide(mol1) and is_alcohol(mol2)) or (is_alkyl_halide(mol2) and is_alcohol(mol1)):
            st.write("- Williamson ether synthesis likely")
        elif is_alkyl_halide(mol1) or is_alkyl_halide(mol2):
            st.write("- Nucleophilic substitution (SN2) possible")

# Template selection dropdown
template_names = list(working_templates.keys())
template_names.insert(0, "Auto-detect")  # Add auto-detect option
selected_template = st.selectbox("Select reaction template", template_names, 
                               help="Choose a specific reaction template or let the system auto-detect")

# Additional option to display only the main product
show_leaving_groups = st.checkbox("Show leaving groups in product", value=False, 
                                 help="When checked, the product will include leaving groups like HBr, H2O, etc.")

# Predict product button
if mol1 and mol2:
    predict_col, reset_col = st.columns([2, 1])
    with predict_col:
        predict_button = st.button("🔬 Predict Product", use_container_width=True)
    with reset_col:
        reset_button = st.button("🔄 Reset All", use_container_width=True)
    
    if predict_button:
        with st.spinner("Predicting product..."):
            try:
                if selected_template == "Auto-detect":
                    # Try all templates and find matches
                    product_dict = predict_reaction_products(mol1, mol2)
                    
                    # Store matches for display
                    st.session_state['template_matches'] = product_dict
                    
                    if product_dict:
                        # Take the first successful template's products
                        first_template = list(product_dict.keys())[0]
                        products = product_dict[first_template]
                        
                        # Join products with dot notation
                        full_product = ".".join(products) if isinstance(products, list) else products
                        reaction_info = f"Matched template: {first_template}"
                    else:
                        st.error("No templates matched these reactants.")
                        full_product = None
                        reaction_info = "No matching templates found"
                else:
                    # Use the specific template selected by the user
                    try:
                        from rxnutils.chem.reaction import ChemicalReaction
                        
                        # Get the selected template
                        template_smarts = working_templates[selected_template]
                        
                        # Create reaction from pattern
                        rxn = ChemicalReaction(template_smarts)
                        rxn.generate_reaction_template(radius=1)
                        
                        # Apply the template to predict products
                        reactants = f"{mol1}.{mol2}"
                        product_list = rxn.canonical_template.apply(reactants)
                        
                        # Flatten the product list
                        products = []
                        if product_list:
                            for product_set in product_list:
                                products.extend(product_set)
                        
                        if products:
                            # Join products with dot notation
                            full_product = ".".join(products)
                            reaction_info = f"Template: {selected_template}"
                            
                            # Store single template match for display
                            st.session_state['template_matches'] = {selected_template: products}
                        else:
                            st.warning(f"The selected template '{selected_template}' didn't match these reactants.")
                            
                            # Try with auto-detect as fallback
                            product_dict = predict_reaction_products(mol1, mol2)
                            
                            if product_dict:
                                st.session_state['template_matches'] = product_dict
                                first_template = list(product_dict.keys())[0]
                                products = product_dict[first_template]
                                full_product = ".".join(products) if isinstance(products, list) else products
                                reaction_info = f"Fell back to template: {first_template}"
                            else:
                                full_product = None
                                reaction_info = "No matching templates found"
                    except Exception as e:
                        st.error(f"Error applying template: {str(e)}")
                        full_product = None
                
                if full_product:
                    # Determine which product to display based on user preference
                    if not show_leaving_groups and "." in full_product:
                        display_product = get_main_product(full_product)
                    else:
                        display_product = full_product
                    
                    # For energy calculations, always use the main product
                    energy_product = get_main_product(full_product) if "." in full_product else full_product
                    
                    # Show the selected product representation
                    st.subheader("Predicted Product")
                    with st.expander("SMILES"):
                        st.code(display_product)
                    st.caption(reaction_info)
                    
                    # Display full reaction if showing leaving groups
                    if show_leaving_groups and "." in full_product:
                        with st.expander("Full Reaction"):
                            st.markdown(f"**Reactants**: {mol1} + {mol2}")
                            st.markdown(f"**Products**: {full_product}")
                            components = full_product.split(".")
                            if len(components) > 1:
                                st.markdown(f"**Main product**: {components[0]}")
                                st.markdown(f"**Leaving group(s)**: {'.'.join(components[1:])}")
                            
                    # Show matching templates for debugging
                    with st.expander("Matching Reaction Templates"):
                        if 'template_matches' in st.session_state and st.session_state['template_matches']:
                            st.write("These templates successfully matched your reactants:")
                            for template, products in st.session_state['template_matches'].items():
                                st.markdown(f"- **{template}**: `{products if isinstance(products, str) else ', '.join(products)}`")
                        else:
                            st.write("No templates matched directly.")
                    
                    # Generate 3D visualization for the display product
                    molblock = generate_3D(display_product)
                    if molblock:
                        with st.expander("3D Product View"):
                            visualize_3D(molblock)
                        
                        try:
                            # Always use the main product for energy calculations
                            E_prod, elements, coords = calculate_energy_with_rdkit(energy_product)
                            xyz_path = "xyz_files/product.xyz"
                            write_xyz_file(elements, coords, xyz_path)
                            with st.expander("Energy"):
                                st.success(f"Product Energy: {E_prod:.6f} Hartree")
                            
                            if "mol1_energy" in st.session_state and "mol2_energy" in st.session_state:
                                result = Energy_comparison(
                                    st.session_state["mol1_energy"], 
                                    st.session_state["mol2_energy"], 
                                    E_prod
                                )
                                
                                st.subheader("Thermodynamic Analysis")
                                E1 = st.session_state["mol1_energy"]
                                E2 = st.session_state["mol2_energy"]
                                delta_E = E1 + E2 - E_prod
                                
                                col1, col2 = st.columns(2)
                                with col1:
                                    st.metric("ΔE (Hartree)", f"{delta_E:.6f}")
                                with col2:
                                    st.metric("ΔE (kcal/mol)", f"{delta_E * 627.5:.2f}")
                                
                                if result == "stable":
                                    st.success("✅ Reaction is thermodynamically favorable at 0 K.")
                                elif result == "equilibrium":
                                    st.info("⚖️ Reaction leads to thermodynamic equilibrium at 0 K.")
                                else:
                                    st.error("❌ Reaction is NOT thermodynamically favorable at 0 K.")
                        except Exception as e:
                            st.error(f"Product energy calculation failed: {e}")
                else:
                    st.error("Could not predict a product for these reactants.")
            except Exception as e:
                st.error(f"Error predicting product: {e}")
                import traceback
                st.code(traceback.format_exc())
    
    if reset_button:
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.rerun()

# Add a section to show available templates
with st.expander("Available Reaction Templates"):
    st.write("These are the reaction templates currently available in the system:")
    for name, smarts in working_templates.items():
        st.markdown(f"**{name}**: `{smarts}`")
    st.info("The app will try to match your reactants against these templates. If none match, it will fall back to more generic reaction prediction methods.")

with st.expander("About this app"):
    st.write("""
    This application predicts chemical reactions and their thermodynamic favorability at 0 K.
    
    **How it works:**
    1. Draw two molecules in the editors above
    2. Select a specific reaction template or let the app auto-detect the reaction
    3. Choose whether to show leaving groups (like HBr, H2O) in the product
    4. Click "Predict Product" to see the reaction product and energy analysis
    5. The app uses RDKit and rxnutils with a template-based approach for chemical predictions
    
    **Features:**
    - Template-based reaction prediction for accurate results
    - Special handling for alkyl halide coupling reactions (Wurtz reaction)
    - Support for specific reaction types like SN2 substitutions
    - Energy calculations using RDKit force fields
    
    **Limitations:**
    - Predictions are based on reaction templates and may not capture all possible reactions
    - Energy calculations are approximate and use 0 K as the reference temperature
    - More complex reactions may not be predicted correctly
    """)