from __future__ import annotations
from typing import Optional
import Bio
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.vectors import Vector, rotmat
import numpy as np


def get_model(filepath: str) -> Optional[Bio.PDB.Model.Model]:
    if filepath.lower().endswith(".cif"):
        return MMCIFParser().get_structure("", filepath)[0]
    elif filepath.lower().endswith(".pdb"):
        return PDBParser().get_structure("", filepath)[0]


def translate_chain(
    model: Bio.PDB.Model.Model, chain_char: str, translate: list[float]
) -> None:
    translation_vec = np.array(translate, "f")

    chain_atoms = model[chain_char].get_atoms()
    for atom in chain_atoms:
        atom.transform(rotmat(Vector(0, 0, 0), Vector(0, 0, 0)), translation_vec)


def rotate_chain(
    model: Bio.PDB.Model.Model, chain_char: str, rotate: list[float]
) -> None:
    rotation_vec = Vector(*rotate)

    chain_atoms = model[chain_char].get_atoms()
    for atom in chain_atoms:
        atom.transform(rotmat(atom.get_vector(), rotation_vec), np.zeros(3))


def save_model(model, filepath) -> None:
    # output in cif file called {original cif name}_transformed.cif
    io = MMCIFIO()
    io.set_structure(model)
    io.save(f'./{filepath.rsplit("/", maxsplit=1)[-1].split(".")[0]}_transformed.cif')
