"""Main logic."""
import sys

from cli import parser
from transformations import transform


def main():
    """Parse args and perform flag processes.

    Get model from args and perform args flag operations on model.

    :return: None
    """
    args = parser.parse_args(sys.argv[1:])

    if model := transform.get_model(args.filepath):
        if args.translate:
            transform.translate_chain(model, args.chain, args.translate)

        if args.rotate:
            transform.rotate_chain(model, args.chain, args.rotate)

        transform.save_model(model, args.filepath)


if __name__ == "__main__":
    main()
