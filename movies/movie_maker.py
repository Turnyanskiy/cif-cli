from pymol import cmd
from glob import glob

lst = glob('./traj_*.cif')
lst.sort()
for traj in lst: cmd.load(traj, 'mov')

print(len(lst))
