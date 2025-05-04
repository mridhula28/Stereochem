import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pub
from rdkit import Chem
from rdkit.Chem import Draw
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from stereochem.generate_isomers import generate_isomers

# ---- Initialize session states ----
if "main_smiles" not in st.session_state:
    st.session_state.main_smiles = ""  # The main molecule input (SMILES)
if "guessed_molecules" not in st.session_state:
    st.session_state.guessed_molecules = set()  # Set to store correct guesses
if "score" not in st.session_state:
    st.session_state.score = 0  # User's score
if "show_answers" not in st.session_state:
    st.session_state.show_answers = False  # Toggle to show answers
if "hint" not in st.session_state:
    st.session_state.hint = False  # Toggle for hint display
if "show_chiral_atoms" not in st.session_state:
    st.session_state.show_chiral_atoms = False

# ---- Update function for the main molecule ----
def update_input_molecule(new_smiles):
    st.session_state.main_smiles = new_smiles
    st.session_state.guessed_molecules = set()  # Reset guesses when new molecule is input
    st.session_state.score = 0
    st.session_state.show_answers = False
    st.session_state.hint = False
    st.session_state.show_chiral_atoms = False

    for key in list(st.session_state.keys()):
        if key.startswith("Atom"):
            st.session_state[key] = False

# ---- Page Config ----
st.set_page_config(page_title="StereoChem", page_icon=":test_tube:", layout="wide")

# ---- Page Title ----
st.title('Stereoisomers in Chemistry')
st.caption("Practical Programming in Chemistry miniproject")

# ---- Tabs for Drawing and Isomers ----
tab1, tab2, tab3 = st.tabs(["Input a molecule", "Draw isomers", "Chirality"])

# ---- Tab 1: Input Molecule (by drawing or name) ----
with tab1:
    st.subheader("Input a Molecule")
    st.markdown("You can input a molecule by drawing it.  \nPlease click on 'Apply' and then on 'Submit Drawing' to input your molecule.")
    
    # --- Drawing input via Ketcher ---
    with st.form(key="draw_form"):
        drawn_smiles = st_ketcher(key="draw_molecule")
        submit_draw = st.form_submit_button("Submit Drawing")

    if submit_draw and drawn_smiles:
        mol_check = Chem.MolFromSmiles(drawn_smiles)
        if mol_check:
            input_2 = Chem.MolToSmiles(mol_check, isomericSmiles=True, canonical=True)
            if input_2 != st.session_state.main_smiles:
                update_input_molecule(input_2)
        else:
            st.warning("Invalid drawn molecule.")

    # --- Sidebar for Molecule Name Input ---
    with st.sidebar.form(key="name_form"):
        molecule_name_input = st.text_input("Enter molecule name:")
        submit_name = st.form_submit_button("Submit Name")

    if submit_name and molecule_name_input:
        try:
            compounds = pub.get_compounds(molecule_name_input, 'name')
            if compounds:
                input_1 = compounds[0].isomeric_smiles
                if input_1 != st.session_state.main_smiles:
                    update_input_molecule(input_1)
            else:
                st.sidebar.warning("Molecule not found.")
        except Exception as e:
            st.sidebar.warning(f"Error fetching molecule: {e}")

    # --- Display the current main molecule ---
    st.sidebar.markdown("### Current Main Molecule:")
    if st.session_state.main_smiles:
        mol = Chem.MolFromSmiles(st.session_state.main_smiles)
        st.sidebar.image(Draw.MolToImage(mol), caption="Main Molecule", use_container_width=True)
        st.sidebar.markdown(f"SMILES: `{st.session_state.main_smiles}`")
    else:
        st.sidebar.info("No molecule selected yet.")

# ---- Tab 2: Show Isomers ----
with tab2:
    st.subheader("Draw and Guess Stereoisomers")
    st.markdown("Draw all possible stereoisomers of the input molecule.  \nPlease click on 'Apply' and then on 'Submit Guess' to submit your drawing.")

    # --- Feedback placeholders ---
    points_placeholder = st.empty()
    message_placeholder = st.empty()

    # --- Generate isomers from main molecule ---
    if st.session_state.main_smiles:
        isomer_set = generate_isomers(st.session_state.main_smiles)

        # --- Form for drawing and guessing ---
        with st.form(key="guess_form"):
            drawn_isomers = st_ketcher(key="draw_isomers")
            submit_isomer = st.form_submit_button("Submit Guess")

            # --- Handle submit guess logic ---
            if submit_isomer and drawn_isomers:
                drawn_mol = Chem.MolFromSmiles(drawn_isomers)
                if drawn_mol:
                    drawn_canon_smiles = Chem.MolToSmiles(drawn_mol, isomericSmiles=True, canonical=True)

                    # Check if drawn molecule is a correct stereoisomer
                    if drawn_canon_smiles in isomer_set and drawn_canon_smiles not in st.session_state.guessed_molecules:
                        st.session_state.guessed_molecules.add(drawn_canon_smiles)
                        st.session_state.score += 1
                        message_placeholder.success("This stereoisomer matches one of the possible stereoisomers!")
                        st.image(Draw.MolToImage(drawn_mol), caption="Drawn Molecule", width=150)
                        st.markdown(f"**Drawn SMILES:** `{drawn_canon_smiles}`")

                        # Check if all stereoisomers have been found
                        if len(st.session_state.guessed_molecules) == len(isomer_set) and len(st.session_state.guessed_molecules) != 0:
                            message_placeholder.success("Congratulations! You found all the stereoisomers!")
                            st.balloons()
                    else:
                        message_placeholder.warning("This stereoisomer is NOT among the generated stereoisomers.")

        # --- Display score ---
        points_placeholder.markdown(f"### Score: {st.session_state.score}")

        # --- Buttons for answer display ---
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button("Show Answers"):
                st.session_state.show_answers = True
        with col2:
            if st.button("Hide Answers"):
                st.session_state.show_answers = False
        with col3:
            if st.button("Hint"):
                st.session_state.hint = True

        # --- Display isomers if answers are shown ---
        if st.session_state.show_answers:
            st.subheader("All Possible Stereoisomers")
            cols = st.columns(4)
            for i, isomer_smiles in enumerate(sorted(isomer_set)):
                mol = Chem.MolFromSmiles(isomer_smiles)
                img = Draw.MolToImage(mol, size=(200, 200))
                col = cols[i % 4]

                try:
                    results = pub.get_compounds(isomer_smiles, namespace='smiles')
                    if results and results[0].iupac_name:
                        iupac_name = results[0].iupac_name
                    else:
                        iupac_name = "Unknown Name"
                except Exception as e:
                    iupac_name = "Unknown Name"

                with col:
                    st.image(img, caption=iupac_name, use_container_width=True)

        # --- Display hint if needed ---
        if st.session_state.hint:
            st.info(f"There are {len(sorted(isomer_set))} possible stereoisomers.")

        # --- Display guessed molecules ---
        if st.session_state.guessed_molecules:
            st.subheader("Your Guessed Stereoisomers")
            guessed_cols = st.columns(4)
            for i, guessed_smiles in enumerate(sorted(st.session_state.guessed_molecules)):
                guessed_mol = Chem.MolFromSmiles(guessed_smiles)
                guessed_img = Draw.MolToImage(guessed_mol, size=(200, 200))
                col = guessed_cols[i % 4]
                with col:
                    st.image(guessed_img, caption=guessed_smiles, use_container_width=True)

    else:
        st.info("Please input a molecule name or draw a molecule first.")

# ---- Tab 3: Chirality ----
with tab3:
    st.subheader("Chirality")
    st.markdown("Find all the chiral centers of the input molecule.")
    message_placeholder = st.empty()
    chiral_atoms = []

    # --- Determine chirality from the molecule ---
    if st.session_state.main_smiles:
        mol = Chem.MolFromSmiles(st.session_state.main_smiles)
        for atom in mol.GetAtoms(): 
            if atom.HasProp('_ChiralityPossible'): 
                chiral_atoms.append(atom.GetIdx())

        img_highlight = Draw.MolToImage(mol, highlightAtoms=chiral_atoms, size=(120,120))

        # --- Display molecule with atom numbers ---
        for atom in mol.GetAtoms():
            atom.SetProp("atomNote", str(atom.GetIdx()))
        img_numbered = Draw.MolToImage(mol, size=(200, 200))
        st.image(img_numbered, caption="Molecule with atom numbers")

        # --- User selection for chiral atoms ---
        number = mol.GetNumAtoms()
        c1, c2 = st.columns(2)
        user_selection = []

        with c1: 
            for i in range(int(number // 2)): 
                if st.checkbox(f"Atom {i}", key=f"Atom_{i}"):
                    user_selection.append(i)

        with c2:
            for i in range(int(number // 2), number):
                if st.checkbox(f"Atom {i}", key=f"Atom_{i}"):
                    user_selection.append(i)


        # --- New: Button for "No chiral atoms" ---
        no_chiral_button = st.button("No chiral atoms")

        # --- Check if user selected correct chiral atoms ---
        if sorted(user_selection) == sorted(chiral_atoms) and chiral_atoms and user_selection:
            message_placeholder.success("Congratulations! You found all the chiral atoms!")
            st.balloons()

        # --- Check case: molecule has no chiral atoms and user clicked the button ---
        elif no_chiral_button and not chiral_atoms:
            message_placeholder.success("Correct! This molecule has no chiral atoms.")
            st.balloons()

        # --- Buttons to show/hide chiral atoms ---
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Show Chiral Atoms"):
                st.session_state.show_chiral_atoms = True
        with col2:
            if st.button("Hide Chiral Atoms"):
                st.session_state.show_chiral_atoms = False

        # --- Display answer if requested ---
        if st.session_state.show_chiral_atoms:
            st.subheader("All chiral atoms highlighted")
            st.image(img_highlight, caption="Chiral Centers Highlighted")
    else:
        st.info("Please input a molecule name or draw a molecule first.")







'''import streamlit as st
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
tab1, tab2, tab3 = st.tabs(["Input a molecule", "Draw isomers", "Chirality"])

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
    st.markdown("""<style> [data-testid="stSidebar"] {background-color: lightblue;} </style>""", unsafe_allow_html=True)
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

    # --- Placeholders for feedback ---
    points_placeholder = st.empty()
    message_placeholder = st.empty()
    subheader_placeholder = st.empty ()
    correct_molecules_placeholder = st.empty()

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
                 # --- IUPAC name ---

                try:
                    results = pub.get_compounds(isomer_smiles, namespace='smiles')
                    if results and results[0].iupac_name:
                        iupac_name = results[0].iupac_name
                    else:
                        iupac_name = "Unknown Name"
                except Exception as e:
                    iupac_name = "Unknown Name"
                with col:
                    st.image(img, caption=iupac_name, use_container_width=True)
       

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
                    
                    message_placeholder.success("This stereoisomer matches one of the possible stereoisomers!")
                    if len(st.session_state.correct_molecules) == len(canon_isomer_set):
                            message_placeholder.success("Congratulations! You found all the stereoisomers!")
                            st.balloons()
                else:
                    message_placeholder.warning("This stereoisomer is NOT among the generated stereoisomers.")
            else:
                message_placeholder.error("Could not parse the drawn SMILES. Please draw a valid molecule.")

        # --- Display current score ---
        points_placeholder.markdown(f"### Score: {st.session_state.score}")
        
        # Display all previously correct molecules horizontally
        if st.session_state.correct_molecules:
            st.subheader("Found Correct Molecules")
    
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

# --- Tab 3: Chirality ---
with tab3: 
    st.subheader("Chirality")
    st.markdown("Find all the chiral centers of the molecule.")
    message_placeholder = st.empty()
    chiral_atoms = []

    # --- Determining chirality from smiles ---
    smiles_input = None
    if st.session_state.input_mode == "name" and st.session_state.molecule_name:
        try:
            compounds = pub.get_compounds(st.session_state.molecule_name, 'name')
            if compounds:
                smiles_input = compounds[0].isomeric_smiles
                mol = Chem.MolFromSmiles(smiles_input)
        except:
            pass
        mol = Chem.MolFromSmiles(smiles_input)
        for atom in mol.GetAtoms(): 
            if atom.HasProp('_ChiralityPossible'): 
                chiral_atoms.append(atom.GetIdx())
        img_highlight = Draw.MolToImage(mol, highlightAtoms=chiral_atoms, size=(120,120))

    elif st.session_state.input_mode == "draw" and "drawn_smiles" in st.session_state:
        smiles_input = st.session_state.drawn_smiles

    # --- Draw the numbered atoms in a molecule ---
    for atom in mol.GetAtoms():
        atom.SetProp("atomNote", str(atom.GetIdx()))
    img = Draw.MolToImage(mol, size=(200, 200))
    st.image(img, caption="Molecule with atom numbers")

    # ---Boxes---
    st.markdown("Tick the boxes whose number corresponds to a chiral atom")
    number = mol.GetNumAtoms()
    c1, c2 = st.columns (2)
    user_selection = []

    with c1: 
        for i in range (0,int(number//2)): 
            if st.checkbox(f"Atom {i}"):
                user_selection.append(i)
    
    with c2:
        for i in range(int(number // 2), number):
            if st.checkbox(f"Atom {i}"):
                user_selection.append(i)

    # --- Check user input ---
    if user_selection == chiral_atoms: 
        message_placeholder.success("Congratulations! You found all the chiral atoms!")
        st.balloons()

    # --- Buttons ---
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Show Chiral Atoms"):
            st.session_state.show_chiral_atoms = True
    with col2:
        if st.button("Hide Chiral Atoms"):
            st.session_state.show_chiral_atoms = False


    if smiles_input:
        # --- Display answer if requested ---
        if st.session_state.show_chiral_atoms:
            st.subheader("All chiral atoms")
            st.image(img_highlight, caption="chiral_atoms")'''
