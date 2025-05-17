from stereochem.functions import generate_isomers

def test_no_stereochemistry():
    assert generate_isomers("CCO") == {"CCO"}

def test_double_bond_stereoisomerism():
    result = generate_isomers("CC=CC")
    assert len(result) == 2

def test_fully_assigned_stereochemistry():
    assert generate_isomers("F[C@H](Cl)Br") == {"F[C@H](Cl)Br"}

def test_chiral_no_isomer():
    # A molecule with a chiral center but too symmetric to produce two distinct stereoisomers
    result = generate_isomers("C[C@@H](C)C")
    assert len(result) == 1

def test_isomer_without_chirality():
    # 2-butene: E/Z isomerism but no chiral center
    result = generate_isomers("CC=CC")
    assert len(result) == 2

def test_no_chirality_no_isomer():
    result = generate_isomers("CCCC")
    assert result == {"CCCC"}

def test_chirality_and_double_bond():
    result = generate_isomers("C/C=C[C@H](Cl)F")
    assert len(result) == 2 or len(result) == 4

def test_invalid_molecule():
    result = generate_isomers("HELP")
    assert result == {"molecule not found"}