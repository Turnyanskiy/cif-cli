from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.vectors import Vector, rotmat
import numpy as np

def transform_chain(filepath, chain, translate, rotate):
    # translation and rotation vectors
    translation = np.array(translate, 'f')
    rotation = Vector(*rotate)

    # getting selected chain
    structure = MMCIFParser().get_structure('', filepath)
    model = structure[0]
    chain = model[chain]
    residues = chain.get_residues()
    
    # performing required transformation on chain
    for residue in residues:
        for atom in residue.get_atoms():
            atom.transform(rotmat(atom.get_vector(), rotation), translation)

    # output in cif file called {orignal cif name}_transformed.cif
    io = MMCIFIO()
    io.set_structure(model)
    io.save(f'./{filepath.rsplit("/", maxsplit=1)[-1].split(".")[0]}_transformed.cif')
