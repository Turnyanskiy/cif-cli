"""Functions used in transformations of pdb/mmCif file."""
from __future__ import annotations

import Bio
import numpy as np
from Bio.PDB.vectors import Vector, rotmat


def translate_chain(
    model: Bio.PDB.Model.Model, chain_char: str, translate: list[float]
) -> None:
    """Translate specified chain from model by vector.

    :param model: The model to transform
    :param chain_char: The model's chain to transform
    :param translate: A vector to translate the chain
    :return: None
    """
    translation_vec = np.array(translate, "f")

    chain_atoms = model[chain_char].get_atoms()
    for atom in chain_atoms:
        atom.transform(rotmat(Vector(0, 0, 0), Vector(0, 0, 0)), translation_vec)


def rotate_chain(
    model: Bio.PDB.Model.Model, chain_char: str, rotate: list[float]
) -> None:
    """Rotate specified chain from model to vector.

    :param model: The model to transform
    :param chain_char: The model's chain to transform
    :param rotate: The vector to rotate the chain onto
    :return: None
    """
    rotation_vec = Vector(*rotate)

    chain_atoms = model[chain_char].get_atoms()
    for atom in chain_atoms:
        atom.transform(rotmat(atom.get_vector(), rotation_vec), np.zeros(3))
