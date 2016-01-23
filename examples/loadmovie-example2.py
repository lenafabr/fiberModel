# load structures as movie into pymol
from pymol import cmd
from glob import *

# get the list of files
files = glob("examples/example2.*.dump.pdb")
# sort the list of files numerically
files.sort( key=lambda x: int(x.split('.')[1]) )

for fname in files:   
    print fname
    if os.path.isfile(fname):
       cmd.load(fname,"struct")
   
   
