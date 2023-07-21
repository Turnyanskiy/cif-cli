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


def save_pose(pose, filepath: str) -> None:
    """Save pose in cif file.

    Pose is saved in current working directory under the name
    of "{old file name}_transformed.cif"

    :param pose: The model to save
    :param filepath: Filepath to original pdb/cif file
    :return: None
    """
    return pose.dump_pdb(
        f'./{filepath.rsplit("/", maxsplit=1)[-1].split(".")[0]}_transformed.cif'
    )
