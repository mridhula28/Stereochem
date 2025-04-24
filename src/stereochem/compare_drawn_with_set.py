# compare the smile of the drawn molecule with the simles of the set of stereoisomers

import streamlit as st
from stereochem.generate_isomers import generate_isomers
from stereochem.example_module import example_module

if smile_code:  # Only proceed if user has drawn something
    drawn_mol = Chem.MolFromSmiles(smile_code)
    
    if drawn_mol:
        # Canonical SMILES with stereochemistry
        drawn_canon_smiles = Chem.MolToSmiles(drawn_mol, isomericSmiles=True, canonical=True)

        # Display drawn structure
        st.image(Draw.MolToImage(drawn_mol), caption="Drawn Molecule")
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
