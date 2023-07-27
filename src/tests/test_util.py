import os

import Bio
import pyrosetta
import pytest

import cif_cli.util.biopython
import cif_cli.util.pyrosetta


@pytest.mark.parametrize(
    "filepath, expected",
    [
        ("./samples/7g93.cif", True),
        ("./samples/7g93.pdb", True),
        ("./samples/7g93", False),
    ],
)
def test_get_model(filepath, expected):
    assert (
        isinstance(cif_cli.util.biopython.get_model(filepath), Bio.PDB.Model.Model)
        == expected
    )


def test_save_model_as_pdb():
    filepath = "samples/7g93.cif"
    model = cif_cli.util.biopython.get_model(filepath)
    cif_cli.util.biopython.save_model_as_pdb(model, "test.pdb")
    os.system("rm -r test.pdb")


def test_save_model_as_cif():
    filepath = "samples/7g93.cif"
    model = cif_cli.util.biopython.get_model(filepath)
    cif_cli.util.biopython.save_model_as_pdb(model, "test.cif")
    os.system("rm -r test.cif")


@pytest.fixture(scope="session")
def test_pyrosetta_init():
    pyrosetta.init(
        "-packing:ex1 "
        "-packing:ex2aro "
        "-run:constant_seed "
        "-score:docking_interface_score true "
        "-flexPepDocking:lowres_preoptimize true "
        "-flexPepDocking:pep_refine true "
        "-flexPepDocking:receptor_chain A "
        "-flexPepDocking:peptide_chain C"
    )


@pytest.mark.usefixtures("test_pyrosetta_init")
@pytest.mark.parametrize(
    "filepath, expected",
    [
        ("./samples/7g93.cif", True),
        ("./samples/7g93.pdb", True),
        ("./samples/7g93", False),
    ],
)
def test_get_pose(filepath, expected):
    assert (
        isinstance(cif_cli.util.pyrosetta.get_pose(filepath), pyrosetta.Pose)
        == expected
    )


@pytest.mark.usefixtures("test_pyrosetta_init")
def test_save_pose_as_pdb():
    filepath = "samples/7g93.cif"
    model = cif_cli.util.pyrosetta.get_pose(filepath)
    cif_cli.util.pyrosetta.save_pose_as_pdb(model, "test.pdb")
    os.system("rm -r test.pdb")
