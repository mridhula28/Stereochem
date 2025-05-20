import unittest
from unittest.mock import patch
import streamlit as st
from stereochem.functions import update_input_molecule


class TestUpdateInputMolecule(unittest.TestCase):

    def setUp(self):
        # Ensure a clean state before each test
        st.session_state.clear()

    def test_empty_state_reset(self):
        update_input_molecule("CCO")
        self.assertEqual(st.session_state.main_smiles, "CCO")
        self.assertEqual(st.session_state.guessed_molecules, set())
        self.assertEqual(st.session_state.score, 0)
        self.assertFalse(st.session_state.show_answers)
        self.assertFalse(st.session_state.hint)
        self.assertFalse(st.session_state.show_chiral_atoms)
        self.assertEqual(st.session_state.chrono_text, "")
        self.assertNotIn("start_time", st.session_state)
        self.assertNotIn("end_time_structures", st.session_state)
        self.assertNotIn("all_iupac_validated", st.session_state)
        self.assertNotIn("balloons_shown", st.session_state)

    def test_partial_tab1_input(self):
        st.session_state.main_smiles = "CCC"
        update_input_molecule("C=CC")
        self.assertEqual(st.session_state.main_smiles, "C=CC")
        self.assertEqual(st.session_state.guessed_molecules, set())

    def test_tab2_filled(self):
        st.session_state.guessed_molecules = {"C[C@H](Cl)F"}
        st.session_state.score = 3
        st.session_state.show_answers = True
        update_input_molecule("CCF")
        self.assertEqual(st.session_state.guessed_molecules, set())
        self.assertEqual(st.session_state.score, 0)
        self.assertFalse(st.session_state.show_answers)

    def test_tab3_filled_checkboxes(self):
        st.session_state["Atom_0"] = True
        st.session_state["Atom_5"] = True
        st.session_state["Atom_9"] = True
        update_input_molecule("CCC")
        self.assertFalse(st.session_state["Atom_0"])
        self.assertFalse(st.session_state["Atom_5"])
        self.assertFalse(st.session_state["Atom_9"])

    def test_all_filled(self):
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

    def test_resets_hint_and_chiral(self):
        st.session_state.update({
            "show_answers": True,
            "hint": True,
            "show_chiral_atoms": True
        })
        update_input_molecule("CCO")
        self.assertFalse(st.session_state.show_answers)
        self.assertFalse(st.session_state.hint)
        self.assertFalse(st.session_state.show_chiral_atoms)

    def test_removes_transient_fields(self):
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

    def test_overwrites_main_smiles(self):
        st.session_state.main_smiles = "old"
        update_input_molecule("new")
        self.assertEqual(st.session_state.main_smiles, "new")


if __name__ == "__main__":
    unittest.main()
