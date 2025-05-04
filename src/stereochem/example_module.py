import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pub
from rdkit import Chem
from rdkit.Chem import Draw
import time
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

    # Reset timing-related variables
    if "end_time_structures" in st.session_state: 
        del st.session_state.end_time_structures
    if "start_time" in st.session_state:
        del st.session_state.start_time
    if "all_iupac_validated" in st.session_state:
        del st.session_state.all_iupac_validated

        # Reset timing and validation states
    for key in ["end_time_structures", "start_time", "all_iupac_validated"]:
        if key in st.session_state:
            del st.session_state[key]

    # Reset dynamic atom checkbox states
    for key in list(st.session_state.keys()):
        if key.startswith("Atom"):
            del st.session_state[key]

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
    score_placeholder = st.empty()
    chrono_placeholder = st.empty()
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
                        import time
                        if "start_time" not in st.session_state:
                            st.session_state.start_time = time.time()
                        message_placeholder.success("This stereoisomer matches one of the possible stereoisomers!")
                        st.image(Draw.MolToImage(drawn_mol), caption="Drawn Molecule", width=150)
                        st.markdown(f"**Drawn SMILES:** `{drawn_canon_smiles}`")

                        # Check if all stereoisomers have been found
                        if len(st.session_state.guessed_molecules) == len(isomer_set) and len(st.session_state.guessed_molecules) != 0:
                            if "end_time_structures" not in st.session_state:
                                import time
                                st.session_state.end_time_structures = time.time()
                                elapsed = st.session_state.end_time_structures - st.session_state.start_time
                                minutes = int(elapsed // 60)
                                seconds = int(elapsed % 60)
                                chrono_placeholder.markdown(f"### Chrono: {minutes} min {seconds} sec")
                                st.balloons()
                            else:
                                message_placeholder.success("🎉 Congratulations! You found all the stereoisomers!")

        # --- Display score ---
        score_placeholder.markdown(f"### Score: {st.session_state.score}")

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
            if "validated_names" not in st.session_state:
                st.session_state.validated_names = set()
            guessed_cols = st.columns(4)
            for i, guessed_smiles in enumerate(sorted(st.session_state.guessed_molecules)):
                guessed_mol = Chem.MolFromSmiles(guessed_smiles)
                guessed_img = Draw.MolToImage(guessed_mol, size=(200, 200))
                col = guessed_cols[i % 4]
                with col:
                    if "name_validation_status" not in st.session_state:
                        st.session_state.name_validation_status = {}

                    with st.form(key=f"form_inline_{guessed_smiles}"):
                        st.image(guessed_img, use_container_width=True)
                        if st.session_state.get(f"reset_requested_{guessed_smiles}", False):
                            st.session_state[f"input_inline_{guessed_smiles}"] = ""
                            st.session_state[f"reset_requested_{guessed_smiles}"] = False
                            st.rerun()
                        user_input_name = st.text_input(
                            f"IUPAC name for isomer {i+1}:", 
                            key=f"input_inline_{guessed_smiles}"
                        )

                        with st.container():
                            st.markdown("<div style='text-align:center;'>", unsafe_allow_html=True)
                            submit_name_check = st.form_submit_button("Check IUPAC Name")
                            reset_name = st.form_submit_button("Reset")
                            st.markdown("</div>", unsafe_allow_html=True)

                        if submit_name_check:
                            try:
                                results = pub.get_compounds(guessed_smiles, namespace='smiles')
                                if results and results[0].iupac_name:
                                    correct_name = results[0].iupac_name.lower().replace("-", "").replace(",", "").replace(" ", "")
                                    user_name = user_input_name.lower().replace("-", "").replace(",", "").replace(" ", "")
                                    if user_name == correct_name:
                                        if guessed_smiles not in st.session_state.validated_names:
                                            st.session_state.score += 1
                                            st.session_state.validated_names.add(guessed_smiles)
                                            score_placeholder.markdown(f"### Score: {st.session_state.score}")
                                        st.session_state.name_validation_status[guessed_smiles] = "correct"
                                    else:
                                        st.session_state.name_validation_status[guessed_smiles] = "incorrect"
                                else:
                                    st.warning("Could not retrieve the correct name from PubChem.")
                            except Exception as e:
                                st.warning(f"An error occurred while retrieving name: {e}")

                        if reset_name:
                            st.session_state[f"reset_requested_{guessed_smiles}"] = True
                            st.session_state.name_validation_status[guessed_smiles] = None
                            st.rerun()

                        if st.session_state.name_validation_status.get(guessed_smiles) == "correct":
                            st.success("Correct IUPAC name!")
                        elif st.session_state.name_validation_status.get(guessed_smiles) == "incorrect":
                            st.error("Incorrect.")

            if "all_iupac_validated" not in st.session_state:
                st.session_state.all_iupac_validated = False

            if (
                not st.session_state.all_iupac_validated
                and st.session_state.guessed_molecules == st.session_state.validated_names
                and len(st.session_state.guessed_molecules) > 0
            ):
                st.balloons()
                st.session_state.all_iupac_validated = True
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

