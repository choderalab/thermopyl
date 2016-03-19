from thermopyl import thermoml_lib
from thermopyl.utils import get_fn
import tempfile
import os, os.path

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

def test_build_pandas_dataframe():
    tmpdir = tempfile.mkdtemp()

    from thermopyl.utils import build_pandas_dataframe, pandas_dataframe

    # Generate dataframe
    filenames = [get_fn("je8006138.xml")]
    [data, compounds] = build_pandas_dataframe(filenames)

    # Write as HDF5
    data.to_hdf(os.path.join(tmpdir, 'data.h5'), 'data')
    compounds.to_hdf(os.path.join(tmpdir, 'compound_name_to_formula.h5'), 'data')

    # Read dataframe
    df = pandas_dataframe(thermoml_path=tmpdir)

    # Clean up tmpdir
    import shutil
    shutil.rmtree(tmpdir)
