from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

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