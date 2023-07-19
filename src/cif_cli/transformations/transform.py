"""Functions used in transformations of pdb/mmCif file."""
from __future__ import annotations

from typing import Optional

import Bio
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.vectors import Vector, rotmat


def get_model(filepath: str) -> Optional[Bio.PDB.Model.Model]:
    """Retrieve the first model from a pdb/cif file.

    :param filepath: Filepath to pdb/cif file
    :return: The first model of the pdb/cif file
    """
    if filepath.lower().endswith(".cif"):
        return MMCIFParser().get_structure("", filepath)[0]
    elif filepath.lower().endswith(".pdb"):
        return PDBParser().get_structure("", filepath)[0]
    else:
        return None


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


def save_model(model, filepath) -> None:
    """Save model in cif file.

    Model is saved in current working directory under the name
    of "{old file name}_transformed.cif"

    :param model: The model to save
    :param filepath: Filepath to original pdb/cif file
    :return: None
    """
    io = MMCIFIO()
    io.set_structure(model)
    io.save(f'./{filepath.rsplit("/", maxsplit=1)[-1].split(".")[0]}_transformed.cif')
