from cli import parser
from transformations import transform

def main():
    args = parser.parse_args()
    transform.transform_chain(args.filepath, args.chain, args.translate, args.rotate)

if __name__ == '__main__':
    main()
