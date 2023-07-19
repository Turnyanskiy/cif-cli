"""Parser for cli inputs."""
from __future__ import annotations

import argparse


def parse_args(args: list[str]) -> argparse.Namespace:
    """Parse args using argparse.

    :param args: Used for tests. Passes flags into function.
    :return: Object containing parsed args.
    """
    # define transform parser
    parser = argparse.ArgumentParser()
    transform_parser = parser.add_subparsers(title="subcommand").add_parser("transform")

    transform_parser.add_argument("filepath")
    transform_parser.add_argument("chain")

    transform_parser.add_argument("-t", "--translate", nargs=3, type=float)
    transform_parser.add_argument("-r", "--rotate", nargs=3, type=float)

    return parser.parse_args(args)
