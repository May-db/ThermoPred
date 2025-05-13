import os
import streamlit as st
from stmol import showmol
from streamlit_ketcher import st_ketcher
import py3Dmol
from pathlib import Path
from reaction_utils import (
    generate_3D, 
    smiles_to_3d, 
    write_xyz_file, 
    calculate_energy_with_rdkit, 
    Energy_comparison, 
    get_product,
    predict_product_with_templates,
    get_main_product,
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
    view.setBackgroundColor('white')
    showmol(view, height=400, width=400)

def draw_and_process(title, session_key):
    """Draw a molecule and calculate its energy."""
    st.subheader(f"{title}")
    smiles = st_ketcher(st.session_state.get(session_key, ""), key=f"{session_key}_ketcher", height=400)
    
    if smiles:
        st.session_state[session_key] = smiles
        st.code(smiles, language="chemical/x-smiles")
        
        molblock = generate_3D(smiles)
        if molblock:
            with st.expander("3D View"):
                visualize_3D(molblock)
            
            try:
                st.info("Calculating energy with RDKit MMFF94/UFF...")
                
                energy, elements, coords = calculate_energy_with_rdkit(smiles)
                
                # register the optimized geometry
                os.makedirs("xyz_files", exist_ok=True)
                xyz_path = os.path.abspath(f"xyz_files/{session_key}_optimized.xyz")
                write_xyz_file(elements, coords, xyz_path)
                
                st.session_state[f"{session_key}_energy"] = energy
                st.success(f"Energy (RDKit): {energy:.6f} Hartree")
                
                # Ajouter une visualisation de la structure optimisée
                #with st.expander("Optimized 3D Structure"):
                    #st.text("Optimized molecular geometry:")
                    # Si vous avez une fonction pour visualiser les fichiers XYZ:
                    # visualize_xyz(xyz_path)
                    # Sinon, afficher le chemin du fichier
                    #st.text(f"Saved to: {xyz_path}")
                
            except Exception as e:
                st.error(f"Energy calculation failed: {str(e)}")
                import traceback
                st.code(traceback.format_exc())
    
    return smiles


if "mol1" not in st.session_state:
    st.session_state.mol1 = ""
if "mol2" not in st.session_state:
    st.session_state.mol2 = ""


col1, col2 = st.columns(2)
with col1:
    mol1 = draw_and_process("Molecule 1", "mol1")
with col2:
    mol2 = draw_and_process("Molecule 2", "mol2")

# Reaction type selection, du coup ca u lieu des datasets
reaction_types = ["Auto-detect", "Substitution", "Addition", "Elimination", "Condensation"]
selected_type = st.selectbox("Select reaction type (optional)", reaction_types)

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
                if selected_type == "Auto-detect":
                    # Get product with or without leaving groups based on user preference
                    full_product = get_product(None, mol1, mol2, return_main_product_only=False)
                else:
                    rxn_type_map = {
                        "Substitution": "substitution",
                        "Addition": "addition",
                        "Elimination": "elimination",
                        "Condensation": "condensation"
                    }
                    rxn_type = rxn_type_map.get(selected_type, "default")
                    
                    # Use the existing predict_product_with_templates function
                    result = predict_product_with_templates(mol1, mol2)
                    if isinstance(result, tuple) and len(result) == 2:
                        full_product, _ = result
                    else:
                        full_product = result
                
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
                    st.code(display_product, language="chemical/x-smiles")
                    
                    # Display full reaction if showing leaving groups
                    if show_leaving_groups and "." in full_product:
                        with st.expander("Full Reaction"):
                            st.markdown(f"**Reactants**: {mol1} + {mol2}")
                            st.markdown(f"**Products**: {full_product}")
                            components = full_product.split(".")
                            if len(components) > 1:
                                st.markdown(f"**Main product**: {components[0]}")
                                st.markdown(f"**Leaving group(s)**: {'.'.join(components[1:])}")
                    
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
    
    if reset_button:
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.rerun()


with st.expander("About this app"):
    st.write("""
    This application predicts chemical reactions and their thermodynamic favorability at 0 K.
    
    **How it works:**
    1. Draw two molecules in the editors above
    2. Select a reaction type (or let the app auto-detect)
    3. Choose whether to show leaving groups (like HBr, H2O) in the product
    4. Click "Predict Product" to see the reaction product and energy analysis
    5. The app uses RDKit and rxnutils for chemical predictions and xTB for energy calculations
    
    **Limitations:**
    - Predictions are based on common reaction patterns and may not capture all possible reactions
    - Energy calculations are approximate and use 0 K as the reference temperature
    - More complex reactions may not be predicted correctly
    """)