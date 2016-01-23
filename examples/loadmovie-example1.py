# load structures as movie into pymol
from pymol import cmd

for f in range(5):   
   fname = 'examples/ex.basinhop.%d.replic.pdb'%(f+1)
   print fname
   if os.path.isfile(fname):
      cmd.load(fname,"struct")
   
   
