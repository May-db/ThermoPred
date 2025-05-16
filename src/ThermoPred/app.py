import os
import streamlit as st
from stmol import showmol
from streamlit_ketcher import st_ketcher
import py3Dmol
from pathlib import Path
from rdkit import Chem

# Import functions from reaction_utils (using only the original product functions)
from reaction_utils import (
    generate_3D, 
    smiles_to_3d, 
    write_xyz_file, 
    calculate_energy_with_rdkit, 
    Energy_comparison, 
    REACTION_TEMPLATES,
    get_main_product,
    react  # Using your original react() function instead of predict_reaction_products
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

st.subheader("Reaction")
# Template selection dropdown
template_names = list(REACTION_TEMPLATES.keys())
template_names.insert(0, "Auto-detect")  # Add auto-detect option
selected_template = st.selectbox("Select reaction template", template_names, 
                               help="Choose a specific reaction template or let the system auto-detect")

# Additional option to display only the main product
#show_leaving_groups = st.checkbox("Show leaving groups in product", value=False, 
                                 #help="When checked, the product will include leaving groups like HBr, H2O, etc.")

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
                    # Use your original react() function for auto-detection
                    product_smiles, template_used = react(mol1, mol2)
                    
                    if product_smiles:
                        reaction_info = f"Matched template: {template_used}"
                        full_product = product_smiles
                    else:
                        st.error("No templates matched these reactants.")
                        full_product = None
                        reaction_info = "No matching templates found"
                else:
                    # Use the specific template selected by the user
                    try:
                        # Use your original react() function with forced template
                        product_smiles, template_used = react(mol1, mol2)
                        
                        if product_smiles and template_used == selected_template:
                            full_product = product_smiles
                            reaction_info = f"Template: {selected_template}"
                        else:
                            st.warning(f"The selected template '{selected_template}' didn't match these reactants.")
                            # Fall back to auto-detect
                            product_smiles, template_used = react(mol1, mol2)
                            if product_smiles:
                                full_product = product_smiles
                                reaction_info = f"Fell back to template: {template_used}"
                            else:
                                full_product = None
                                reaction_info = "No matching templates found"
                    except Exception as e:
                        st.error(f"Error applying template: {str(e)}")
                        full_product = None
                
                if full_product:
                    # Determine which product to display based on user preference
                    #if not show_leaving_groups and "." in full_product:
                    #if  "." in full_product:
                        #display_product = get_main_product(full_product)
                    #else:
                    display_product = full_product
                    
                    # For energy calculations, always use the main product
                    #energy_product = get_main_product(full_product) if "." in full_product else full_product
                    energy_product = full_product
                    
                    # Show the selected product representation
                    st.subheader("Predicted Product")
                    st.caption(reaction_info)
                    productketcher = st_ketcher(full_product, key="product_ketcher", height=400)
                    with st.expander("SMILES"):
                        st.code(display_product)
                    
                    # Display full reaction if showing leaving groups
                    #if show_leaving_groups and "." in full_product:
                        #with st.expander("Full Reaction"):
                            #st.markdown(f"**Reactants**: {mol1} + {mol2}")
                            #st.markdown(f"**Products**: {full_product}")
                            #components = full_product.split(".")
                            #if len(components) > 1:
                                #st.markdown(f"**Main product**: {components[0]}")
                                #st.markdown(f"**Leaving group(s)**: {'.'.join(components[1:])}")
                    
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
    for name, smarts in REACTION_TEMPLATES.items():
        st.markdown(f"**{name}**: `{smarts}`")
    st.info("The app will try to match your reactants against these templates.")

with st.expander("About this app"):
    st.write("""
    This application predicts chemical reactions and their thermodynamic favorability at 0 K.
    
    **Features:**
    - Uses template-based reaction prediction from reaction_utils.py
    - Energy calculations using RDKit force fields
    - 3D visualization of reactants and products
    
    **Note:** This version uses only the original product prediction functions from reaction_utils.py
    """)