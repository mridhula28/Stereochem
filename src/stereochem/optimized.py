# Refactored version of the Streamlit-based stereochemistry tool
# with clear modular functions for better structure and reuse

import os
import sys
import base64
import time
import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pub
from rdkit import Chem
from rdkit.Chem import Draw

# Local module imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from stereochem.functions import generate_isomers, update_input_molecule

# ------------------------- Configuration and Setup -------------------------

def set_page():
    st.set_page_config(page_title="StereoChem", page_icon=":test_tube:", layout="wide")

def set_background(image_path):
    def encode_image(path):
        with open(path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode()
    encoded_img = encode_image(image_path)
    page_bg_img = f"""
    <style>
    .stApp {{
        background-image: url("data:image/png;base64,{encoded_img}");
        background-size: cover;
        background-repeat: no-repeat;
        background-attachment: fixed;
    }}
    </style>
    """
    st.markdown(page_bg_img, unsafe_allow_html=True)

def initialize_session():
    defaults = {
        "main_smiles": "",
        "guessed_molecules": set(),
        "score": 0,
        "show_answers": False,
        "hint": False,
        "show_chiral_atoms": False,
        "chrono_text": ""
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

# ------------------------- Molecule Input Tab -------------------------

def draw_input_section():
    st.subheader("Input a Molecule")
    st.markdown("Draw a molecule. Click 'Apply' and then 'Submit Drawing'.")
    with st.form(key="draw_form"):
        drawn_smiles = st_ketcher(key="draw_molecule")
        submitted = st.form_submit_button("Submit Drawing")
    if submitted and drawn_smiles:
        mol = Chem.MolFromSmiles(drawn_smiles)
        if mol:
            Chem.RemoveStereochemistry(mol)
            smiles = Chem.MolToSmiles(mol, canonical=True)
            if smiles != st.session_state.main_smiles:
                update_input_molecule(smiles)
        else:
            st.warning("Invalid drawn molecule.")

def name_input_section():
    with st.sidebar.form(key="name_form"):
        name = st.text_input("Enter molecule name:")
        submit = st.form_submit_button("Submit Name")
    if submit and name:
        try:
            results = pub.get_compounds(name, 'name')
            if results:
                mol = Chem.MolFromSmiles(results[0].isomeric_smiles)
                if mol:
                    Chem.RemoveStereochemistry(mol)
                    smiles = Chem.MolToSmiles(mol, canonical=True)
                    if smiles != st.session_state.main_smiles:
                        update_input_molecule(smiles)
            else:
                st.sidebar.warning("Molecule not found.")
        except Exception as e:
            st.sidebar.warning(f"Error fetching molecule: {e}")

def display_main_molecule():
    st.sidebar.markdown("### Current Main Molecule:")
    smiles = st.session_state.main_smiles
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        st.sidebar.image(Draw.MolToImage(mol), caption="Main Molecule", use_container_width=True)
        st.sidebar.markdown(f"SMILES: `{smiles}`")
    else:
        st.sidebar.info("No molecule selected yet.")
    if st.sidebar.button("🔄 Reset Everything"):
        update_input_molecule("")

# ------------------------- Draw Isomers Tab -------------------------

def draw_isomers_tab():
    st.subheader("Draw Isomers")
    if not st.session_state.main_smiles:
        st.info("Please select or draw a molecule first.")
        return

    isomer_smiles_list = generate_isomers(st.session_state.main_smiles)
    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("### Guess the isomers")
        with st.form("guess_form"):
            user_guess = st_ketcher(key="isomer_guess")
            guess_submitted = st.form_submit_button("Submit Guess")

        if guess_submitted and user_guess:
            try:
                mol = Chem.MolFromSmiles(user_guess)
                if mol:
                    guessed = Chem.MolToSmiles(Chem.RemoveStereochemistry(Chem.Mol(mol)), canonical=True)
                    if guessed in isomer_smiles_list and guessed not in st.session_state.guessed_molecules:
                        st.session_state.guessed_molecules.add(guessed)
                        st.session_state.score += 1
                        st.success("Correct isomer!")
                    elif guessed in st.session_state.guessed_molecules:
                        st.warning("Already guessed.")
                    else:
                        st.error("Not a correct isomer.")
            except:
                st.error("Invalid SMILES string.")

    with col2:
        st.metric("Score", st.session_state.score)
        if st.checkbox("Show answers"):
            st.session_state.show_answers = True

    if st.session_state.show_answers:
        st.markdown("### All Stereoisomers")
        for smi in isomer_smiles_list:
            mol = Chem.MolFromSmiles(smi)
            st.image(Draw.MolToImage(mol), caption=smi, use_container_width=True)

# ------------------------- Chirality Tab -------------------------

def chirality_tab():
    st.subheader("Chirality Check")
    smiles = st.session_state.main_smiles
    if not smiles:
        st.info("Please select or draw a molecule first.")
        return

    mol = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    st.markdown("### Molecule Structure with Chiral Centers")
    st.image(Draw.MolToImage(mol, highlightAtoms=[a.GetIdx() for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]), use_container_width=True)

    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if chiral_centers:
        st.markdown("**Chiral Centers:**")
        for idx, chirality in chiral_centers:
            st.markdown(f"- Atom index {idx}: {chirality}")
    else:
        st.info("No chiral centers detected.")

# ------------------------- Page Setup -------------------------

def main():
    set_page()
    set_background("assets/background.png")
    initialize_session()

    st.title("Stereoisomers in Chemistry")
    st.caption("Practical Programming in Chemistry miniproject")

    tab1, tab2, tab3 = st.tabs(["Input a molecule", "Draw isomers", "Chirality"])

    with tab1:
        draw_input_section()
        name_input_section()
        display_main_molecule()

    with tab2:
        draw_isomers_tab()

    with tab3:
        chirality_tab()

if __name__ == "__main__":
    main()
