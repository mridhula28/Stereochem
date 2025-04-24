# test/manual_test_generate_isomers.py

import sys
import os

# Add the 'src' folder to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

# Import the function from stereochem in src
from stereochem.generate_isomers import generate_isomers

def test_generate_isomers():
    smiles = "CC(O)C(=O)O"  # Example: lactic acid
    
    # Call the generate_isomers function
    isomers = generate_isomers(smiles)
    
    # Print the result
    print("Generated isomers:", isomers)

# Run the test manually
if __name__ == "__main__":
    test_generate_isomers()
