import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pub
from rdkit import Chem 
from rdkit.Chem import Draw

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from stereochem.generate_isomers import generate_isomers

#setting page title and icon
st.set_page_config(page_title= "StereoChem", page_icon= ":test_tube:", layout= "wide") 

#Head setter 
st.title('Stereoisomers in Chemistry')
st.caption("Practical Proramming In Chemistry miniproject")
st.markdown("Draw all possible stereoisomers of the inputed molecule")


#molecule drawing frame 
molecule_name = st.sidebar.text_input("Enter the name of the molecule")
if molecule_name:
    try:
        compounds = pub.get_compounds(molecule_name, 'name')
        if compounds:
            smiles_mol = compounds[0].isomeric_smiles
            isomer_set_RS_EZ = generate_isomers(smiles_mol)
            # Display the result (afficher les smiles de tous les isomeres du set)
            # st.sidebar.markdown("### Generated Stereoisomers")
            # for i, isomer in enumerate(isomer_set_RS_EZ, 1):
            #     st.sidebar.text(f"{i}: {isomer}")
            mol = Chem.MolFromSmiles(smiles_mol)
            st.sidebar.image(Draw.MolToImage(mol), caption=f"{molecule_name} structure")
            st.sidebar.markdown(f"Smile code: {smiles_mol}")
        else:
            st.warning("No compound found. Please check the molecule name.")
    except Exception as e:
        st.error(f"Error fetching data: {e}")
smile_code = st_ketcher(molecule_name)


