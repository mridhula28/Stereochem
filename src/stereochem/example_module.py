import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pub
from rdkit import Chem 
from rdkit.Chem import Draw
from typing import Any
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from stereochem.generate_isomers import generate_isomers

score = 0
#setting page title and icon
st.set_page_config(page_title= "StereoChem", page_icon= ":test_tube:", layout= "wide") 

isomer_set_RS_EZ = None

#Head setter, placeholders
st.title('Stereoisomers in Chemistry')
st.caption("Practical Programming In Chemistry miniproject")
points = st.empty() # placeholder for point 
image = st.empty() # placeholder for images 
message = st.empty() # placeholder for messages 
st.markdown("Draw all possible stereoisomers of the inputed molecule")

def scores (score, points): 
    points.markdown(f"Score = {score}")
    return score


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

#Display answers
if "show_answers" not in st.session_state:
    st.session_state["show_answers"] = False

#Display hints 
if "Hint" not in st.session_state: 
    st.session_state["Hint"] = False 

#Create two columns to place the two buttons side by side
col1, col2, col3, col4 = st.columns(4)

with col1:
    if st.button("Show Answers"):
        st.session_state["show_answers"] = True

with col2:
    if st.button("Hide Answers"):
        st.session_state["show_answers"] = False

with col3: 
    if st.button ("Hint"): 
        st.session_state["Hint"] = True    

with col4: 
    if st.button("Hide hint"): 
        st.session_state["Hint"] = False 


#Display isomers if enabled
if molecule_name and isomer_set_RS_EZ and st.session_state["show_answers"]:
    st.subheader("All possible stereoisomers")
    cols = st.columns(4)
    for i, isomer_smiles in enumerate(sorted(isomer_set_RS_EZ)):
        mol = Chem.MolFromSmiles(isomer_smiles)
        img = Draw.MolToImage(mol, size=(200, 200))
        col = cols[i % 4]
        with col:
            st.image(img, caption=isomer_smiles, use_container_width=True)

# gives a hint one the number of stereoisomers there needs to be found 
if molecule_name and isomer_set_RS_EZ and st.session_state["Hint"]: 
    st.subheader(f"There is a total of {len(sorted(isomer_set_RS_EZ))} stereoisomers")


if smile_code:  # Only proceed if user has drawn something
    drawn_mol = Chem.MolFromSmiles(smile_code)
    
    if drawn_mol:
        # Canonical SMILES with stereochemistry
        drawn_canon_smiles = Chem.MolToSmiles(drawn_mol, isomericSmiles=True, canonical=True)

        # Compare with generated stereoisomers
        if molecule_name and 'isomer_set_RS_EZ' in locals():
            # Normalize all generated isomers as canonical SMILES
            canon_isomer_set = {Chem.MolToSmiles(Chem.MolFromSmiles(sm), isomericSmiles=True, canonical=True) for sm in isomer_set_RS_EZ}
            
            if drawn_canon_smiles in canon_isomer_set:
                score = score + 1
                print(score)
                message.success("This stereoisomer matches one of the possible stereoisomers.")
                image.image(Draw.MolToImage(drawn_mol), caption="Drawn Molecule", width=100)
                st.markdown(f"**Drawn SMILES:** `{drawn_canon_smiles}`")
                scores(score, points)
                
            else:
                message.warning("This stereoisomer is NOT among the generated stereoisomers.")
        else:
            message.info("Stereoisomer set not ready yet.")
    else:
        message.error("Could not parse the drawn SMILES. Please draw a valid molecule.")



