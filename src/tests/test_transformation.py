import copy

from Bio.PDB.vectors import Vector, rotmat

import cif_cli.util.biopython
from cif_cli.transformations import transform


def test_translate_chain():
    original_model = cif_cli.util.biopython.get_model("samples/7g93.cif")
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
    original_model = cif_cli.util.biopython.get_model("samples/7g93.cif")
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
