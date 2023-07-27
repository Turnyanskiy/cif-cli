"""Module containing implementation of FlexPepDock docking protocol.

This protocol is intended for docking between flexible peptides and globular proteins.

When using FlexPepDockProtocol, please ensure Rosetta has been initialised with the
appropriate options - see the below link. If the correct options are not initialised,
FlexPepDockProtocol may not impact pose conformation.

FlexPepDock docking is composed of the following steps:
    - Initial sidechain packing
    - FlexPepDock peptide to MHC
    - Final sidechain packing

For information on the FlexPepDock protocol, see:
https://www.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock

Typical usage example:

    import pyrosetta

    from tcr.protocols.pmhc_gen import FlexPepDockProtocol

    # Initialise Rosetta
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

    # Load a pose into memory. Peptide chain id is C, while MHC id is A
    pose = pyrosetta.pose_from_pdb(...)

    # Dock it using FlexPepDockProtocol
    FlexPepDockProtocol(docking_partners="A_C").apply(pose=pose)
"""
import pyrosetta
from pyrosetta import Vector1
from pyrosetta.rosetta.core import pack
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover


class FlexPepDockProtocol:
    """Implementation of FlexPepDock (FPD) flexible peptide docking protocol.

    This protocol is intended for docking between flexible peptides and globular
    proteins.

    When using FlexPepDockProtocol, please ensure Rosetta has been initialised with the
    appropriate options - see the below link. If the correct options are not initialised,
    FlexPepDockProtocol may not impact pose conformation.

    For information on the FlexPepDock protocol, see:
    https://www.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock

    Attributes:
        docking_partners: String with format 'A..C_X..Z' representing <protein_peptide> chains
            to FlexPepDock dock. Chains either side of the '_' will be docked together.
    """

    def __init__(self, docking_partners: str) -> None:
        """Instantiate a FlexPepDockProtocol object.

        Args:
            docking_partners: String with format 'A..C_X..Z' representing
            <protein_peptide> chains to FlexPepDock dock. Chains either side of the '_'
            will be docked together.
        """
        self.docking_partners = docking_partners

    def __repr__(self) -> str:
        """Return a string representation of a FlexPepDockProtocol instance."""
        docking_partners = self.docking_partners
        return f"{type(self).__name__}({docking_partners})"

    def apply(self, pose: Pose) -> None:
        """Apply FlexPepDock protocol a PyRosetta pMHC pose. This function has side-effects.

        FlexPepDock docking is composed of the following steps:
            - Initial sidechain packing
            - FlexPepDock peptide to MHC
            - Final sidechain packing

        Args:
            pose: PyRosetta pMHC pose.
        """
        setup_foldtree(pose, self.docking_partners, Vector1([1]))
        self._pack_sidechains(pose=pose)
        self._flex_pep_dock(pose=pose)
        self._pack_sidechains(pose=pose)

    def _pack_sidechains(self, pose: Pose) -> None:
        """Apply PyRosetta sidechain packing to pose. This function has-side effects.

        Args:
            pose: PyRosetta pose.
        """
        pack_score_fn = ScoreFunctionFactory.create_score_function("docking")
        new_traj_energy = self._enable_trajectory(pack_score_fn)
        packer = build_packer(pack_score_fn)
        packer.apply(pose)


    def _enable_trajectory(self, score_function):
        dt_st = pyrosetta.rosetta.core.scoring.ScoreType.dump_trajectory
        score_function.set_weight(dt_st, 1)
        emo = pyrosetta.rosetta.core.scoring.methods.EnergyMethodOptions()
        emo.dump_trajectory_prefix('./test/work')
        emo.dump_trajectory_stride(1000)
        DumpTrajectoryEnergy = pyrosetta.rosetta.core.energy_methods.DumpTrajectoryEnergy
        new_traj_energy = DumpTrajectoryEnergy(emo)
        score_function.all_methods().append(new_traj_energy)
        return new_traj_energy


    def _flex_pep_dock(self, pose: Pose) -> None:
        """Apply PyRosetta FlexPepDock docking to pose. This function has side-effects.

        Args:
            pose: PyRosetta pose.
        """
        flexpepdock_protocol = FlexPepDockingProtocol(1)
        flexpepdock_protocol.apply(pose)


def build_packer(score_function) -> PackRotamersMover:
    """Generate PyRosetta PackRotamersMover for sidechain packing.

    Args:
        score_function: Scoring function for packer.

    Returns:
        Packer object.
    """
    task_factory = pack.task.TaskFactory()
    task_factory.push_back(pack.task.operation.InitializeFromCommandline())
    task_factory.push_back(pack.task.operation.RestrictToRepacking())
    packer = PackRotamersMover(score_function)
    packer.task_factory(task_factory)
    return packer
