from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO


parser = MMCIFParser()

models = parser.get_structure("phmc", "traj_20230724132202274_pack.clean.pdb").get_models()

print(models)


for i, model in enumerate(models):
    print(model, models[i])
    io=MMCIFIO()
    io.set_structure(models[i])    
    io.save(f"{i}.cif")
