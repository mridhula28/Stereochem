# test/manual_test_generate_isomers.py

# Import the function from stereochem in src
from src.stereochem.generate_isomers import generate_isomers

def test_generate_isomers():
    smiles = "CC(O)C(=O)O"  # Example: lactic acid
    
    # Call the generate_isomers function
    isomers = generate_isomers(smiles)
    
    # Print the result
    print("Generated isomers:", isomers)

# Run the test manually
if __name__ == "__main__":
    test_generate_isomers()
