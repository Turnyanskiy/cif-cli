"""Util functions for BioPython."""
from typing import Optional

import Bio
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.mmcifio import MMCIFIO


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


def save_model(model: Bio.PDB.Model.Model, filepath: str) -> None:
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
