from thermopyl import thermoml_lib
from thermopyl.utils import get_fn

formula = "C3H5N2OClBr"
reference_atom_count = 13
reference_element_counts = dict(C=3, H=5, N=2, O=1, Cl=1, Br=1)


def test_count():
    n = thermoml_lib.count_atoms(formula)
    assert n == reference_atom_count


def test_count_atoms_in_set():
    n = thermoml_lib.count_atoms_in_set(formula, ["C", "H"])
    assert n == 8


def test_formula_to_element_counts():
    element_counts = thermoml_lib.formula_to_element_counts(formula)
    assert element_counts == reference_element_counts


def test_thermopyl():
    filename = get_fn("je8006138.xml")

    parser = thermoml_lib.Parser(filename)
    current_data = parser.parse()
    name_to_formula = parser.compound_name_to_formula
