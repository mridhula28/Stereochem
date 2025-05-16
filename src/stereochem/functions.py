from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import streamlit as st

# creates a set containing all the isomers for a given molecule; both R/S and E/Z.

def generate_isomers(smiles: str) -> set:
    # Convert the SMILES string into an RDKit molecule object
    mol_gen_iso = Chem.MolFromSmiles(smiles)
    
    if mol_gen_iso is None:
        return {"molecule not found"}
    
    # Set up options for stereoisomer enumeration:
    # - onlyUnassigned=True: only generate isomers for undefined stereocenters
    # - unique=True: ensures only unique stereoisomers are returned (avoids duplicates)
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True)
    
    # Generate a list of all possible stereoisomers using the specified options
    isomers = list(EnumerateStereoisomers(mol_gen_iso, options=opts))
    
    # Convert each stereoisomer molecule back into a canonical SMILES string with stereochemistry
    # Store them in a set to automatically remove duplicates
    return {Chem.MolToSmiles(iso, isomericSmiles=True) for iso in isomers}

# ---- Update function for the main molecule ----

def update_input_molecule(new_smiles):
    st.session_state.main_smiles = new_smiles
    st.session_state.guessed_molecules = set()
    st.session_state.score = 0
    st.session_state.show_answers = False
    st.session_state.hint = False
    st.session_state.show_chiral_atoms = False
    st.session_state.validated_names = set()
    st.session_state.name_validation_status = {}
    st.session_state.all_iupac_validated = False
    st.session_state.balloons_shown = False
    st.session_state.start_time = None
    st.session_state.end_time_structures = None
    st.session_state.chrono_text = ""

    # Reset chrono
    for key in ["start_time", "end_time_structures"]:
        if key in st.session_state:
            del st.session_state[key]

    # Reset validation
    for key in ["validated_names", "all_iupac_validated", "balloons_shown"]:
        if key in st.session_state:
            del st.session_state[key]

    # Reset atom selection checkboxes
    for key in list(st.session_state.keys()):
        if key.startswith("Atom"):
            st.session_state[key] = False