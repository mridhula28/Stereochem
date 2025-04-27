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

if "correct_molecules" not in st.session_state:
    st.session_state.correct_molecules = []

# --- Page Title ---
st.title('Stereoisomers in Chemistry')
st.caption("Practical Programming In Chemistry miniproject")


# --- Tabs for Drawing and Isomers ---
tab1, tab2 = st.tabs(["Input a molecule", "Draw isomers"])

# --- Tab 1: Input by Drawing ---
with tab1:
    st.subheader("Draw a Molecule")
    drawn_smiles = st_ketcher(key="draw_molecule_tab1")

    if drawn_smiles:
        st.session_state.drawn_smiles = drawn_smiles
        st.session_state.input_mode = "draw"
        st.session_state.score = 0 # reset score to 0 if a new molecule is drawn
        st.session_state.molecule_name = ""  # <- clear the NAME properly!!
        st.session_state.correct_molecules = []

# --- Sidebar for Molecule Input ---
with st.sidebar:
    st.title("Molecule Input")

    if "molecule_name" not in st.session_state:
        st.session_state.molecule_name = ""

    molecule_name_input = st.text_input("Enter the name of the molecule:", value=st.session_state.molecule_name)

    if molecule_name_input != st.session_state.molecule_name:
        st.session_state.molecule_name = molecule_name_input
        st.session_state.score = 0 # reset score to 0 if a new molecule is drawn
        if molecule_name_input:
            st.session_state.input_mode = "name"
            st.session_state.drawn_smiles = ""  # clear drawn smiles if user typed a name
            st.session_state.correct_molecules = []

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



# --- Tab 2: Show Isomers ---
with tab2:
    st.subheader("Isomers")
    st.markdown("Draw all possible stereoisomers of the input molecule.")

    # --- Drawing input ---
    drawn_isomers = st_ketcher(key="draw_molecule_tab2")

    # --- Determine SMILES to use ---
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

    # --- Initialize session states ---
    if "show_answers" not in st.session_state:
        st.session_state.show_answers = False

    if "hint" not in st.session_state:
        st.session_state.hint = False

    if "score" not in st.session_state:
        st.session_state.score = 0

    # --- Placeholders for feedback ---
    points_placeholder = st.empty()
    message_placeholder = st.empty()
    image_placeholder = st.empty()

    # --- Buttons ---
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        if st.button("Show Answers"):
            st.session_state.show_answers = True
    with col2:
        if st.button("Hide Answers"):
            st.session_state.show_answers = False
    with col3:
        if st.button("Hint"):
            st.session_state.hint = True
    with col4:
        if st.button("Hide Hint"):
            st.session_state.hint = False

    # --- Main logic ---
    if smiles_input:
        isomer_set_RS_EZ = generate_isomers(smiles_input)

        # --- Display stereoisomers if requested ---
        if st.session_state.show_answers:
            st.subheader("All Possible Stereoisomers")
            cols = st.columns(4)
            for i, isomer_smiles in enumerate(sorted(isomer_set_RS_EZ)):
                mol = Chem.MolFromSmiles(isomer_smiles)
                img = Draw.MolToImage(mol, size=(200, 200))
                col = cols[i % 4]
                with col:
                    st.image(img, caption=isomer_smiles, use_container_width=True)

        # --- Display hint if requested ---
        if st.session_state.hint:
            st.info(f"There are {len(sorted(isomer_set_RS_EZ))} possible stereoisomers.")

        # --- Handle user drawn molecule checking ---
        if drawn_isomers:
            drawn_mol = Chem.MolFromSmiles(drawn_isomers)

            if drawn_mol:
                # Canonical SMILES of drawn molecule
                drawn_canon_smiles = Chem.MolToSmiles(drawn_mol, isomericSmiles=True, canonical=True)

                # Normalize generated isomers
                canon_isomer_set = {Chem.MolToSmiles(Chem.MolFromSmiles(sm), isomericSmiles=True, canonical=True)
                                    for sm in isomer_set_RS_EZ}

                # Check if drawn structure matches any isomer
                if drawn_canon_smiles in canon_isomer_set:
                    if drawn_canon_smiles not in st.session_state.correct_molecules:
                        st.session_state.correct_molecules.append(drawn_canon_smiles)
                        st.session_state.score = len(st.session_state.correct_molecules)  # "st.session_state.score += 1" was not working
                        
                        if len(st.session_state.correct_molecules) == len(canon_isomer_set):
                            st.success("Congratulations! You found all the stereoisomers! 🎉")
                            st.balloons()
                    
                    message_placeholder.success("This stereoisomer matches one of the possible stereoisomers!")
                    image_placeholder.image(Draw.MolToImage(drawn_mol), caption="Drawn Molecule", width=150)
                    st.markdown(f"**Drawn SMILES:** `{drawn_canon_smiles}`")
                else:
                    message_placeholder.warning("This stereoisomer is NOT among the generated stereoisomers.")
            else:
                message_placeholder.error("Could not parse the drawn SMILES. Please draw a valid molecule.")

        # --- Display current score ---
        points_placeholder.markdown(f"### Score: {st.session_state.score}")
        
        # Display all previously correct molecules horizontally
        if st.session_state.correct_molecules:
            st.subheader("Previously Correct Molecules")
    
            # Create columns dynamically based on the number of correct molecules
            num_molecules = len(st.session_state.correct_molecules)

            if num_molecules > 0:
                cols = st.columns(num_molecules)  # Create as many columns as there are molecules

                for i, correct_smiles in enumerate(st.session_state.correct_molecules):
                    mol = Chem.MolFromSmiles(correct_smiles)
                    img = Draw.MolToImage(mol, size=(120, 120))  # You can adjust the size

                # Display the molecule in the corresponding column
                    with cols[i]:
                        st.image(img, caption=correct_smiles, width=120)

    else:
        st.info("Please input a molecule name or draw a molecule first.")
