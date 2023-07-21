"""Main logic."""
import sys

import pyrosetta
import util.biopython
import util.pyrosetta
from cli import parser
from flexpepdock import flexpepdock
from transformations import transform


def main():
    """Parse args and perform flag processes.

    Get model from args and perform args flag operations on model.
    If flexpepdock is called PyRosetta is initialised.

    :return: None
    """
    args = parser.parse_args(sys.argv[1:])

    if args.subcommand == "flexpepdock":
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

        if pose := util.pyrosetta.get_pose(args.filepath):
            flexpepdock.FlexPepDockProtocol(args.docking_partners).apply(pose)
            util.pyrosetta.save_pose(pose, args.filepath)

    elif args.subcommand == "transform":
        if model := util.biopython.get_model(args.filepath):
            if args.translate:
                transform.translate_chain(model, args.chain, args.translate)

            if args.rotate:
                transform.rotate_chain(model, args.chain, args.rotate)

            util.biopython.save_model(model, args.filepath)


if __name__ == "__main__":
    main()
