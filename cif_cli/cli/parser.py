import argparse

def parse_args():
    # define transform parser
    parser = argparse.ArgumentParser()
    transform_parser = parser.add_subparsers(title='subcommand').add_parser('transform')

    
    transform_parser.add_argument('filepath')
    transform_parser.add_argument('chain')

    transform_parser.add_argument('-t', '--translate', default=[0, 0, 0], nargs=3, type=float)
    transform_parser.add_argument('-r', '--rotate', default=[0, 0, 0], nargs=3, type=float)
    
    return parser.parse_args()


