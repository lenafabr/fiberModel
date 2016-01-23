#!/usr/bin/python

import sys, os, shutil
from DNAparse import *

if len(sys.argv) < 3:
    print """This script grabs the geometric parameters from a pdb file and writes them out to a param file. If the param file already has SPIRALBEND, FILLBEAD, or  POSM, XAXM, etc. lines, then the old ones will be removed. Original parameter file gets copied into paramfile.prev.
    usage: geomFromPDB.py pdbfile paramfile <-u unwrap> <-b DNAres> <-o outfile>
    -u specifies an integer number of basepairs to unwrap from each end when getting the end positions and tangents
    -b specifies which residues in the pdb correspond as DNA; by default, it looks for DA,DT,DC,DG. List the residue names separated by commas (no spaces).    
    -o specifies an output pdb file to check how well the geometric parameters match. Default is pdbfile-mod.pdb"""
    
    sys.exit()

pdbfile = sys.argv[1]
paramfile = sys.argv[2]

if '-u' in sys.argv:
    ind = sys.argv.index('-u')
    offset = int(sys.argv[ind+1])
else:
    offset = 0

if '-b' in sys.argv:
    ind = sys.argv.index('-b')
    DNAres = sys.argv[ind+1].strip().split(',')
else:
    DNAres = ['DA','DT','DG','DC']

if '-o' in sys.argv:
    ind = sys.argv.index('-o')
    outfile = sys.argv[ind+1]
else:
    outfile = os.path.split(os.path.splitext(pdbfile)[0])[-1] + '-mod.pdb'

# set up the full structure and the DNA structure
struct = Structure()
struct.structFromFile(pdbfile)
DNA = DNAobj()
DNA.DNAFromStruct(struct,DNAres)

# PCA on the DNA coordinates
(eval,evec) = DNA.PCA()

# rotate entire structure to line up z axis with the smallest eigenval
struct.rotateM(linalg.inv(evec))
struct.recenter()

axatoms = DNA.getAxis() # central axis of DNA
[comatm,posm,posp] = DNA.getEndPos(array([0,0,0]),offset) # edge positions

# rotate to make edge average on the x axis
zax = array([0,0,1])
wantx = (array(posp)+array(posm))/2;
wantx = wantx/sqrt(dot(wantx,wantx))
wantx = wantx - dot(wantx,zax)*zax
curx = array([1,0,0])
cang = dot(wantx,curx)/sqrt(dot(wantx,wantx))
sgn = dot(cross(curx,wantx),zax)
sang = -sign(sgn)*sqrt(1-cang**2)
rotmat = array([[cang,-sang,0],[sang,cang,0],[0,0,1]])
struct.rotateM(rotmat)

maxz = DNA.atoms[0].coords[2]; minz = DNA.atoms[0].coords[2]

for c in DNA.chains:
    for r in c.residues[4:-4]:
        for a in r.atoms:
            if a.coords[2] > maxz:
                maxz = a.coords[2]
            elif a.coords[2] < minz:
                minz = a.coords[2]                            

axatoms = DNA.getAxis() # central axis of DNA
[tanatoms,tanvecs] = DNA.getTangents(offset,bpback=5) #edge tangents
[comatm,posm,posp] = DNA.getEndPos(array([0,0,0]),offset) # edge positions
# convert positions to nm
posm = 0.1*posm; posp = 0.1*posp

height = maxz-minz
radius = mean([sqrt(a[0]**2 + a[1]**2) for a in DNA.axis[3:-3]])
print 'DNA height, center radius:', height, radius


# get length of DNA
contlen = DNA.getLength(offset)

# output the geometry struff into a pdb file
DNA.atoms.extend(axatoms)
DNA.atoms.extend(tanatoms)
DNA.atoms.append(comatm)
# axes for the structure
#Zatm = Atom(coords=com+[0,0,10],name='COM')
#Xatm = Atom(coords=com+[10,0,0],name='COM')
#Yatm = Atom(coords=com+[0,10,0],name='COM')
#DNA.atoms.extend([Xatm,Yatm,Zatm])

DNA.atoms.append(Atom(coords=[0,0,0],name='AC'))
DNA.atoms.append(Atom(coords=[0,0,-height/2],name='AC1'))
DNA.atoms.append(Atom(coords=[0,0,height/2],name='AC2'))

DNA.renumAtoms()
DNA.outputPDB(outfile)

struct.atoms.extend(axatoms)
struct.atoms.extend(tanatoms)
struct.atoms.append(comatm)

struct.atoms.append(Atom(coords=[0,0,0],name='AC'))
struct.atoms.append(Atom(coords=[0,0,-height/2],name='AC1'))
struct.atoms.append(Atom(coords=[0,0,height/2],name='AC2'))
struct.renumAtoms()
struct.outputPDB(outfile)


# output geometric parameters into param file
shutil.copyfile(paramfile,paramfile+'.prev')
PF = open(paramfile)
lines = PF.readlines()
PF.close()
PF = open(paramfile,'w')
for l in lines:
    if l[0] != '#' and ('SPIRALBENDS' in l or 'BENDTANM' in l \
                        or 'BENDTANP' in l  or 'POSM' in l or 'POSP' in l \
                        or 'XAXM' in l or 'XAXP' in l or 'FILLBEAD' in l):
        pass    
    else:
        PF.write(l.strip()+'\n')

PF.write('POSM %20.10fD0 %20.10fD0 %20.10fD0\n' %(posm[0],posm[1],posm[2]))
PF.write('POSP %20.10fD0 %20.10fD0 %20.10fD0\n' %(posp[0],posp[1],posp[2]))
PF.write('BENDTANM %20.10fD0 %20.10fD0 %20.10fD0\n' %(tanvecs[0][0],tanvecs[0][1],tanvecs[0][2]))
PF.write('BENDTANP %20.10fD0 %20.10fD0 %20.10fD0\n' %(tanvecs[2][0],tanvecs[2][1],tanvecs[2][2]))
PF.write('XAXM %20.10fD0 %20.10fD0 %20.10fD0\n' %(tanvecs[1][0],tanvecs[1][1],tanvecs[1][2]))
PF.write('XAXP %20.10fD0 %20.10fD0 %20.10fD0\n' %(tanvecs[3][0],tanvecs[3][1],tanvecs[3][2]))

# write in all the filler beads
for i in range(offset,len(DNA.axis)-offset):
    bead = (DNA.axis[i])*0.1
    branchatm = DNA.bpairs[i].residues[0].getAtomByName('^C1',1)
    branch = (branchatm.coords)*0.1
    PF.write('FILLBEAD %8.3fD0%8.3fD0%8.3fD0%8.3fD0%8.3fD0%8.3fD0\n'\
             %(bead[0], bead[1],bead[2],branch[0],branch[1],branch[2]))
    
PF.close()


