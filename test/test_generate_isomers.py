import unittest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from stereochem.functions import generate_isomers

class TestGenerateIsomers(unittest.TestCase):

    def test_generate_isomers(self):
        # No stereochemistry
        self.assertEqual(generate_isomers("CCO"), {"CCO"})

        # Double bond: E/Z isomerism
        result = generate_isomers("CC=CC")
        self.assertEqual(len(result), 2)

        # Only geometric isomerism (no chirality)
        result = generate_isomers("CC=CC")
        self.assertEqual(len(result), 2)

        # No chirality or isomerism
        self.assertEqual(generate_isomers("CCCC"), {"CCCC"})

        # Two chiral centers
        result = generate_isomers("FC(Br)C(Cl)I")
        self.assertEqual(len(result), 4)

        # Chiral center + double bond
        result = generate_isomers("CC=CC(Cl)F")
        self.assertEqual(len(result), 4)

        # Invalid molecule
        self.assertEqual(generate_isomers("HELP"), {"molecule not found"})

if __name__ == '__main__':
    unittest.main()
