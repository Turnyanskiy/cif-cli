# `cif-cli`

## Installation

```
./install.sh
```

## Usage

```
Usage:
  python3 src/cif_cli [<path>] [subcommand]

Subcommand:
  flexpepdock <docking_partners>     dock peptide
  transform <chain>                  transform chain

Arguments:
  <path>                 The path to the cif/pdb file
  <chain>                The chain to transform 
  <docking_partners>     The docking partners written in the form <chain>_<chain>

Options:
  -h, --help                            show this help message and exit

  transform:
    -t, --translate FLOAT FLOAT FLOAT     translates the chain by vector
    -r, --rotate FLOAT FLOAT FLOAAT       rotates the chain to the vector
```
