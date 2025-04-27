import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pub
from rdkit import Chem 
from rdkit.Chem import Draw

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from stereochem.generate_isomers import generate_isomers

# --- Page Config ---
st.set_page_config(page_title="StereoChem", page_icon=":test_tube:", layout="wide")

# --- Initialize session states ---
if "input_mode" not in st.session_state:
    st.session_state.input_mode = "name"

if "show_answers" not in st.session_state:
    st.session_state.show_answers = False

# --- Page Title ---
st.title('Stereoisomers in Chemistry')
st.caption("Practical Programming In Chemistry miniproject")

# --- Sidebar for Molecule Input ---
with st.sidebar:
    st.title("Molecule Input")

    if "molecule_name" not in st.session_state:
        st.session_state.molecule_name = ""

    molecule_name_input = st.text_input("Enter the name of the molecule:", value=st.session_state.molecule_name)

    if molecule_name_input != st.session_state.molecule_name:
        st.session_state.molecule_name = molecule_name_input
        if molecule_name_input:
            st.session_state.input_mode = "name"
            st.session_state.drawn_smiles = ""  # clear drawn smiles if user typed a name


    # --- Display the molecule preview dynamically based on input_mode ---
    if st.session_state.input_mode == "name" and st.session_state.molecule_name:
        try:
            compounds = pub.get_compounds(st.session_state.molecule_name, 'name')
            if compounds:
                smiles_mol = compounds[0].isomeric_smiles
                mol = Chem.MolFromSmiles(smiles_mol)
                st.image(Draw.MolToImage(mol), caption=f"{st.session_state.molecule_name} structure")
                st.markdown(f"SMILES: {smiles_mol}")
            else:
                st.warning("No compound found. Please check the molecule name.")
        except Exception as e:
            st.error(f"Error fetching molecule: {e}")

    elif st.session_state.input_mode == "draw" and "drawn_smiles" in st.session_state and st.session_state.drawn_smiles:
        try:
            drawn_mol = Chem.MolFromSmiles(st.session_state.drawn_smiles)
            st.image(Draw.MolToImage(drawn_mol), caption="Drawn Molecule Structure")
            st.markdown(f"SMILES: {st.session_state.drawn_smiles}")
        except Exception as e:
            st.error(f"Error displaying drawn molecule: {e}")


# --- Tabs for Drawing and Isomers ---
tab1, tab2 = st.tabs(["Input a molecule", "Draw isomers"])

# --- Tab 1: Input by Drawing ---
with tab1:
    st.subheader("Draw a Molecule")
    drawn_smiles = st_ketcher(key="draw_molecule_tab1")

    if drawn_smiles:
        st.session_state.drawn_smiles = drawn_smiles
        st.session_state.input_mode = "draw"
        st.session_state.molecule_name = ""  # <- clear the NAME properly!!


# --- Tab 2: Show Isomers ---
with tab2:
    st.subheader("Isomers")
    st.markdown("Draw all possible stereoisomers of the input molecule.")

    drawn_isomers = st_ketcher(key="draw_molecule_tab2")

    # Determine which SMILES to use
    smiles_input = None
    if st.session_state.input_mode == "name" and st.session_state.molecule_name:

        try:
            compounds = pub.get_compounds(st.session_state.molecule_name, 'name')
            if compounds:
                smiles_input = compounds[0].isomeric_smiles
        except:
            pass
    elif st.session_state.input_mode == "draw" and "drawn_smiles" in st.session_state:
        smiles_input = st.session_state.drawn_smiles

    if smiles_input:
        isomer_set_RS_EZ = generate_isomers(smiles_input)

        col1, col2 = st.columns(2)
        with col1:
            if st.button("Show Answers"):
                st.session_state.show_answers = True
        with col2:
            if st.button("Hide Answers"):
                st.session_state.show_answers = False

        if st.session_state.show_answers:
            st.subheader("All Possible Stereoisomers")
            cols = st.columns(4)
            for i, isomer_smiles in enumerate(sorted(isomer_set_RS_EZ)):
                mol = Chem.MolFromSmiles(isomer_smiles)
                img = Draw.MolToImage(mol, size=(200, 200))
                col = cols[i % 4]
                with col:
                    st.image(img, caption=isomer_smiles, use_container_width=True)
    else:
        st.info("Please input a molecule name or draw a molecule first.")