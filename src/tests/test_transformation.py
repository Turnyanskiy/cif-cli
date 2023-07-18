import os
import pytest
from cif_cli.transformations import transform
import Bio
from Bio.PDB.vectors import Vector, rotmat
import copy


@pytest.mark.parametrize(
    "filepath, expected",
    [
        ("./test_transformation/7g93.cif", True),
        ("./test_transformation/7g93.pdb", True),
        ("./test_transformation/7g93", False),
    ],
)
def test_get_model(filepath, expected):
    assert isinstance(transform.get_model(filepath), Bio.PDB.Model.Model) == expected


def test_translate_chain():
    original_model = transform.get_model("./test_transformation/7g93.cif")
    new_model = copy.deepcopy(original_model)
    transform.translate_chain(new_model, "A", [1, 1, 1])

    original_chain_atoms = original_model["A"].get_atoms()
    new_chain_atoms = new_model["A"].get_atoms()

    success = True
    for (original, new) in zip(original_chain_atoms, new_chain_atoms):
        if any(original.get_vector().get_array()) != any(
            (new.get_vector() - Vector(1, 1, 1)).get_array()
        ):
            success = False
            break

    assert success


def test_rotate_chain():
    original_model = transform.get_model("./test_transformation/7g93.cif")
    new_model = copy.deepcopy(original_model)
    transform.rotate_chain(new_model, "A", [1, 1, 1])

    original_chain_atoms = original_model["A"].get_atoms()
    new_chain_atoms = new_model["A"].get_atoms()

    success = True
    for (original, new) in zip(original_chain_atoms, new_chain_atoms):
        rotation_vec = rotmat(new.get_vector(), Vector(1, 1, 1))

        if any(original.get_vector().get_array()) != any(
            (new.get_vector().left_multiply(rotation_vec)).get_array()
        ):
            success = False
            break

    assert success


def test_save_model():
    filepath = "./test_transformation/7g93.cif"
    model = transform.get_model(filepath)
    transform.save_model(model, filepath)
    os.system("rm -r 7g93_transformed.cif")
