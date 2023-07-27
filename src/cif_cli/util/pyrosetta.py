"""Util functions for PyRosetta."""
from typing import Optional

import pyrosetta


def get_pose(filepath: str) -> Optional[pyrosetta.Pose]:
    """Retrieve pose from a pdb/cif file.

    :param filepath: Filepath to pdb/cif file
    :return: The first model of the pdb/cif file
    """
    try:
        print(pyrosetta.pose_from_file(filepath))
        return pyrosetta.pose_from_file(filepath)
    except Exception as e:
        print(e)

    return None


def save_pose_as_pdb(pose, output: str) -> None:
    """Save pose in pdb file.

    Pose is saved using output filepath

    :param pose: The model to save
    :param output: Filepath to output pdb file
    :return: None
    """
    return pose.dump_pdb(output)
