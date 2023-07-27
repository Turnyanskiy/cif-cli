"""Util functions for BioPython."""
from typing import Optional

import Bio
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import PDBIO

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


def save_model_as_cif(model: Bio.PDB.Model.Model, output: str) -> None:
    """Save model in cif file.

    Model is saved using output filepath.

    :param model: The model to save
    :param output: Filepath to output cif file
    :return: None
    """
    io = MMCIFIO()
    io.set_structure(model)
    io.save(output) #f'./{filepath.rsplit("/", maxsplit=1)[-1].split(".")[0]}_transformed.cif'


def save_model_as_pdb(model: Bio.PDB.Model.Model,  output: str) -> None:
    """Save model in pdb file.

    Model is saved using output filepath.

    :param model: The model to save
    :param output: Filepath to output pdb file
    :return: None
    """
    io = PDBIO()
    io.set_structure(model)
    io.save(output)
