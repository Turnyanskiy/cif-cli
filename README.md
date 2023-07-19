# `cif-cli`

## Installation

```
./install.sh
```

## Usage

```
Usage:
  python cif_cli transform [options] [<path>] [<chain>]

Arguments:
  <path> The path to the cif/pdb file for which the transformation should take place on
  <chain> The chain letter for which the transofmration should take place on

Options:
  -h, --help                            show this help message and exit
  -t, --translate FLOAT FLOAT FLOAT     translates the chain by vector
  -r, --rotate FLOAT FLOAT FLOAAT       rotates the chain to the vector
```
