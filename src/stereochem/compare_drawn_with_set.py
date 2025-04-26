# compare the smile of the drawn molecule with the simles of the set of stereoisomers

import streamlit as st
from rdkit import Chem 
from rdkit.Chem import Draw
from streamlit_ketcher import st_ketcher
from stereochem.generate_isomers import generate_isomers
from stereochem.example_module import example_module

# compare the smile of the drawn molecule with the simles of the set of stereoisomers

smile_code = st_ketcher(molecule_name)

if smile_code:                                              # It only proceeds if user has drawn something
    st.session_state["drawn_smiles"] = smile_code

if "drawn_smiles" in st.session_state:
    drawn_mol = Chem.MolFromSmiles(st.session_state["drawn_smiles"])
    
    if drawn_mol:                                                           # Canonical SMILES with stereochemistry
        drawn_canon_smiles = Chem.MolToSmiles(drawn_mol, isomericSmiles=True, canonical=True)
        st.image(Draw.MolToImage(drawn_mol), caption="Drawn Molecule")      # Display drawn structure
        st.markdown(f"**Drawn SMILES:** `{drawn_canon_smiles}`")

        # Compare with generated stereoisomers
        if molecule_name and 'isomer_set_RS_EZ' in locals():
            if drawn_canon_smiles in isomer_set_RS_EZ:
                st.success("This stereoisomer matches one of the possible stereoisomers.")
            else:
                st.warning("This stereoisomer is NOT among the generated stereoisomers.")
        else:
            st.info("Stereoisomer set not ready yet.")
    else:
        st.error("Could not parse the drawn SMILES. Please draw a valid molecule.")
