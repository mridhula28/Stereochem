# compare the smile of the drawn molecule with the simles of the set of stereoisomers

import streamlit as st
from rdkit import Chem 
from rdkit.Chem import Draw
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
        if molecule_name and isomer_set_RS_EZ:
            canon_isomer_set = [
                (Chem.MolFromSmiles(smile), Chem.MolToSmiles(Chem.MolFromSmiles(smile), isomericSmiles=True, canonical=True))
                for smile in isomer_set_RS_EZ
            ]

            mols = []
            highlight = []
            for mol, canon_smile in canon_isomer_set:
                if drawn_canon_smiles == canon_smile:
                    mols.append(mol)
                    highlight.append(True)  # Correct match
                else:
                    mols.append(mol)
                    highlight.append(False)  # No match

            # Now draw them with highlight
            legends = ["Matched" if h else "Not Matched" for h in highlight]
            mols_per_row = 4  # How many molecules in a row

            # Generate a nice grid image
            img = Draw.MolsToGridImage(
                mols,
                molsPerRow=mols_per_row,
                subImgSize=(200, 200),
                legends=legends,
                highlightColor=(0.0, 1.0, 0.0),  # Green
                highlightAtomLists=[list(range(mol.GetNumAtoms())) if matched else [] for mol, matched in zip(mols, highlight)]
            )

            st.image(img, caption="Generated Isomers (Matched Highlighted)")
            else:
                st.warning("This stereoisomer is NOT among the generated stereoisomers.")
        else:
            st.info("Stereoisomer set not ready yet.")
    else:
        st.error("Could not parse the drawn SMILES. Please draw a valid molecule.")