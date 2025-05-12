import streamlit as st
import pubchempy as pub
from rdkit import Chem
from rdkit.Chem import Draw

DEFAULT_SESSION_STATE = {
    "main_smiles": "",
    "guessed_molecules": set(),
    "score": 0,
    "show_answers": False,
    "hint": False,
    "show_chiral_atoms": False,
    "balloons_shown": False,
    "validated_names": set(),
    "name_validation_status": {},
    "all_iupac_validated": False,
}

def reset_states(defaults):
    for key, value in defaults.items():
        st.session_state[key] = value

    for key in list(st.session_state.keys()):
        if key.startswith("Atom_"):
            st.session_state[key] = False


def initialize_session():
    for key, default_value in DEFAULT_SESSION_STATE.items():
        if key not in st.session_state:
            st.session_state[key] = default_value

    for key in list(st.session_state.keys()):
        if key.startswith("Atom_"):
            st.session_state[key] = False

def reset_main_molecule():
    reset_states(DEFAULT_SESSION_STATE)
    st.rerun()

def display_in_grid(items, render_fn, columns=4, **kwargs):
    cols = st.columns(columns)
    for i, item in enumerate(items):
        with cols[i % columns]:
            render_fn(item, i, **kwargs)


def render_guess(guessed_smiles, index, score_placeholder):

    guessed_mol = Chem.MolFromSmiles(guessed_smiles)
    guessed_img = Draw.MolToImage(guessed_mol, size=(200, 200))

    if "name_validation_status" not in st.session_state:
        st.session_state.name_validation_status = {}

    with st.form(key=f"form_inline_{guessed_smiles}"):
        st.image(guessed_img, use_container_width=True)

        if st.session_state.get(f"reset_requested_{guessed_smiles}", False):
            st.session_state[f"input_inline_{guessed_smiles}"] = ""
            st.session_state[f"reset_requested_{guessed_smiles}"] = False
            st.rerun()

        user_input_name = st.text_input(
            f"IUPAC name for isomer {index+1}:", 
            key=f"input_inline_{guessed_smiles}"
        )

        submit_name_check = st.form_submit_button("Check IUPAC Name")
        reset_name = st.form_submit_button("Reset")

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
            reset_name_validation(guessed_smiles)

        status = st.session_state.name_validation_status.get(guessed_smiles)
        if status == "correct":
            st.success("Correct IUPAC name!")
        elif status == "incorrect":
            st.error("Incorrect.")

def reset_name_validation(smiles):
    st.session_state[f"input_inline_{smiles}"] = ""
    st.session_state[f"reset_requested_{smiles}"] = False
    st.session_state.name_validation_status[smiles] = None
    st.rerun()

def get_iupac_name(smiles):
    try:
        results = pub.get_compounds(smiles, namespace='smiles')
        return results[0].iupac_name if results and results[0].iupac_name else "Unknown Name"
    except Exception:
        return "Unknown Name"