import unittest
import streamlit as st
from stereochem.functions import update_input_molecule

class TestUpdateInputMolecule(unittest.TestCase):

    def setUp(self):
        st.session_state.clear()

    def test_update_input_molecule_all_cases(self):
        # --- Test 1: Reset state when starting from empty ---
        update_input_molecule("CCO")
        self.assertEqual(st.session_state.main_smiles, "CCO")
        self.assertEqual(st.session_state.guessed_molecules, set())
        self.assertEqual(st.session_state.score, 0)
        self.assertFalse(st.session_state.show_answers)        # reset to False - meaning no answers shown when a new molecule is input
        self.assertFalse(st.session_state.hint)                # reset to False - meaning no hints when a new molecule is input
        self.assertFalse(st.session_state.show_chiral_atoms)   # reset to False - meaning no chiral atoms highlighted when a new molecule is input
        self.assertEqual(st.session_state.chrono_text, "")     # here the timer display is cleared
        self.assertNotIn("start_time", st.session_state)
        self.assertNotIn("end_time_structures", st.session_state)
        self.assertNotIn("all_iupac_validated", st.session_state)
        self.assertNotIn("balloons_shown", st.session_state)

        # --- Test 2: Partial input update resets guessed molecules ---
        st.session_state.main_smiles = "CCC"
        update_input_molecule("C=CC")
        self.assertEqual(st.session_state.main_smiles, "C=CC")
        self.assertEqual(st.session_state.guessed_molecules, set())

        # --- Test 3: Reset guessed molecules, score, show_answers when Tab2 filled ---
        st.session_state.guessed_molecules = {"C[C@H](Cl)F"}
        st.session_state.score = 3
        st.session_state.show_answers = True
        update_input_molecule("CCF")
        self.assertEqual(st.session_state.guessed_molecules, set())
        self.assertEqual(st.session_state.score, 0)
        self.assertFalse(st.session_state.show_answers)

        # --- Test 4: Reset atom checkboxes ---
        st.session_state["Atom_0"] = True
        st.session_state["Atom_5"] = True
        st.session_state["Atom_9"] = True
        update_input_molecule("CCC")
        self.assertFalse(st.session_state["Atom_0"])
        self.assertFalse(st.session_state["Atom_5"])
        self.assertFalse(st.session_state["Atom_9"])

        # --- Test 5: Full reset (after that all the tabs were filled)  ---
        st.session_state.update({
            "main_smiles": "CCC",
            "guessed_molecules": {"C[C@H](Br)F"},
            "score": 5,
            "show_answers": True,
            "hint": True,
            "show_chiral_atoms": True,
            "validated_names": {"C[C@H](Br)F"},
            "name_validation_status": {"C[C@H](Br)F": "correct"},
            "all_iupac_validated": True,
            "balloons_shown": True,
            "start_time": 123456.0,
            "end_time_structures": 123459.0,
            "chrono_text": "Elapsed",
            "Atom_0": True,
            "Atom_1": True
        })
        update_input_molecule("CCF")
        self.assertEqual(st.session_state.main_smiles, "CCF")
        self.assertEqual(st.session_state.guessed_molecules, set())
        self.assertEqual(st.session_state.score, 0)
        self.assertFalse(st.session_state["Atom_0"])
        self.assertFalse(st.session_state["Atom_1"])
        self.assertEqual(st.session_state.chrono_text, "")
        self.assertNotIn("start_time", st.session_state)
        self.assertNotIn("end_time_structures", st.session_state)
        self.assertNotIn("validated_names", st.session_state)
        self.assertNotIn("all_iupac_validated", st.session_state)
        self.assertNotIn("balloons_shown", st.session_state)

        # --- Test 6: Reset hint and chiral atom flags ---
        st.session_state.update({
            "show_answers": True,
            "hint": True,
            "show_chiral_atoms": True
        })
        update_input_molecule("CCO")
        self.assertFalse(st.session_state.show_answers)
        self.assertFalse(st.session_state.hint)
        self.assertFalse(st.session_state.show_chiral_atoms)

        # --- Test 7: Ensure timing and validation fields are cleared on input update ---
        st.session_state.update({
            "start_time": 123,
            "end_time_structures": 456,
            "validated_names": {"valid"},
            "all_iupac_validated": True,
            "balloons_shown": True
        })
        update_input_molecule("CCC")
        for key in ["start_time", "end_time_structures", "validated_names", "all_iupac_validated", "balloons_shown"]:
            self.assertNotIn(key, st.session_state)

        # --- Test 8: Overwrite main_smiles with new input ---
        st.session_state.main_smiles = "CCN"
        update_input_molecule("CCO")
        self.assertEqual(st.session_state.main_smiles, "CCO")

         # --- Test 9: Re-entering the same molecule should still reset the state ---

        # Set up session state with existing molecule and some filled fields
        st.session_state.update({
            "main_smiles": "CCO",
            "guessed_molecules": {"C[C@H](Br)F"},
            "score": 4,
            "show_answers": True,
            "hint": True,
            "show_chiral_atoms": True,
            "chrono_text": "Elapsed",
            "Atom_0": True,
            "start_time": 123,
            "end_time_structures": 456,
            "validated_names": {"C[C@H](Br)F"},
            "all_iupac_validated": True,
            "balloons_shown": True
        })
    
        # Re-enter the same molecule
        update_input_molecule("CCO")
    
        # Even though it is the same molecule, state should be reset
        self.assertEqual(st.session_state.main_smiles, "CCO")
        self.assertEqual(st.session_state.guessed_molecules, set())
        self.assertEqual(st.session_state.score, 0)
        self.assertFalse(st.session_state.show_answers)
        self.assertFalse(st.session_state.hint)
        self.assertFalse(st.session_state.show_chiral_atoms)
        self.assertEqual(st.session_state.chrono_text, "")
        self.assertFalse(st.session_state["Atom_0"])
    
        # Transient fields should be removed
        for key in ["start_time", "end_time_structures", "validated_names", "all_iupac_validated", "balloons_shown"]:
            self.assertNotIn(key, st.session_state)


if __name__ == "__main__":
    unittest.main()

