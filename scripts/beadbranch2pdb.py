#!/usr/bin/python
# convert a structure file containing list of beads and branches
# as output by fibermodel.exe into a pdb with a helical DNA representation
# for visualization purposes.

from pdbutils import *
import sys, os
from genutils import *
from copy import copy

class Bend:
    # bend object with tails and charges
    def __init__(self,cent,orient,tails = None,charges=None,markpoints=None):
        # initiate object with a given center and orientation
        # tails is a list of lists of coordinates triples
        if tails == None:
            self.tails = []
        else:
            self.tails = tails

        if charges == None:
            self.charges = []
        else:
            self.charges = charges

        if markpoints == None:
            self.markpoints =[]
        else:
            self.markpoints = markpoints

        self.cent = cent        
        self.orient = orient
        

def placeBP(chain,branch,atmtypes,nbp=1,fillnbp=None,twist=0,dnawidth = 2,helang=140,connectmax=-1):
    """ place nbp base pairs between every pair of beads in the chain
     use branch to specify orientation
     returns list of atoms
     dnawidth is the width of the dsDNA in nm
     helang is the angle between corresponding atoms on the two helices
     (measured in degrees around the helix axis).
     The twist is used to figure out in which direction to interpolate beads
     Connect filler specifies whether the filler edges and the linker DNA should be connected."""  

    atoms1 = []; atoms2 = []; atoms3 = []; axes=[]; count = 0
    bpcount = 0;
    for b in range(0,len(chain)-1):
        diff = [chain[b+1][i]-chain[b][i] for i in range(3)]
        ndiff = norm(diff)
        switchdir = (ndiff*twist > pi)
        switchdir = 0;        

        if (atmtypes[b] == 'A' and atmtypes[b+1] == 'A'):
            nbase = nbp
        else:
            nbase = fillnbp
        
        isfiller = (atmtypes[b] == 'F')
        
        #do not place the last linker bead since it is redundant
        if (not isfiller and b < len(atmtypes) and atmtypes[b+1]=='F'):
            continue

        bpspace = 1.0/nbase
        
        # get the tangent and normalize
        tvec = diff
        tvec = normalize(tvec)        
        # get the x axis from the branches
        # x axes at the beads
        xax1 = [branch[b][i] - chain[b][i] for i in range(3)]
        if b < len(chain)-2:
            xax2 = [branch[b+1][i] - chain[b+1][i] for i in range(3)]
        else:
            xax2[:] = xax1[:]
        yax1 = crossprod(tvec,xax1)
        xax1 = normalize(xax1); xax2 = normalize(xax2); yax1 = normalize(yax1)
        phi12 = dihedral(branch[b], chain[b], chain[b+1], branch[b+1])

        # align angle to be close to what we expect it to be
        if phi12 - ndiff*twist > pi:
            phi12 = phi12-2*pi
        elif phi12-ndiff*twist <-pi:
            phi12 = phi12 + 2*pi                       

        if switchdir:
            phi12 = 2*pi+phi12
        
        for bp in range(nbase):
            bpcount = bpcount + 1
            # get the center of the basepair
            if isfiller:
                center = [chain[b][i] for i in range(3)]
            else:
                center = [chain[b][i]+diff[i]*(bpspace/2+bp*bpspace) for i in range(3)]                    
            #center = [chain[b][i]+diff[i]*(bp*bpspace) for i in range(3)]    
                    
            count += 1
            if isfiller:
                name = 'F3'
            else:
                name = 'A3'

            #atoms3.append(Atom(coords=[chain[b][i]+diff[i]*(bp*bpspace) for i in range(3)],num=count,name=name))
            atoms3.append(Atom(coords=center,num=count,name=name))
            if len(atoms3) > 1:                
                if (connectmax<0 or norm(array(atoms3[-1].coords)-array(atoms3[-2].coords)) < connectmax):
                    atoms3[-1].conect.append(atoms3[-2])
                    atoms3[-2].conect.append(atoms3[-1])
            
            # x axis at this particular basepair
            # need to interpolate in angle between bead x axes
            if (b+2 < len(atmtypes) and atmtypes[b] == 'A' and atmtypes[b+1] == 'A' and atmtypes[b+2] == 'F'):
                phi = (bp-(nbase-1)/2.0)*phi12/((nbp+fillnbp)/2)
            elif isfiller:
                phi = 0
            else:
                phi = (bp-(nbase-1)/2.0)*phi12/nbase #+ phi12/nbase/2

            xax = [cos(phi)*xax1[i] + sin(phi)*yax1[i] for i in range(3)]
            yax = crossprod(tvec,xax)
            yax = normalize(yax)
            
            # place two atoms in the x-y plane (first one on the x-axis)  separated by angle helang       
            count += 1
            theta = helang*pi/180              
            coords = [center[i] + xax[i] for i in range(3)]
            if isfiller:
                name = 'F1'
            else:
                name = 'A1'
            atoms1.append(Atom(coords=coords,num=count,name=name))

            count += 1            
            coords = [center[i] + dnawidth/2*(cos(theta)*xax[i] + sin(theta)*yax[i]) for i in range(3)]
            if isfiller:
                name = 'F2'
            else:
                name = 'A2'
            atoms2.append(Atom(coords=coords,num=count,name=name))

            atoms1[-1].resnum = bpcount
            atoms2[-1].resnum = bpcount
            atoms3[-1].resnum = bpcount          

            # keep track of atom connections
            #atoms1[-1].conect.append(atoms2[-1])
            #atoms2[-1].conect.append(atoms1[-1])
            atoms1[-1].conect.append(atoms3[-1])
            atoms2[-1].conect.append(atoms3[-1])
            atoms3[-1].conect.extend([atoms1[-1],atoms2[-1]])

            if len(atoms1) > 1:                
                if (connectmax<0 or norm(array(atoms1[-1].coords)-array(atoms1[-2].coords)) < connectmax):
                    atoms1[-1].conect.append(atoms1[-2])
                    atoms1[-2].conect.append(atoms1[-1])
            if len(atoms2) > 1:
                if (connectmax<0 or norm(array(atoms2[-1].coords)-array(atoms2[-2].coords)) < connectmax):
                    atoms2[-1].conect.append(atoms2[-2])
                    atoms2[-2].conect.append(atoms2[-1])

            
    return atoms1+atoms2+atoms3

def placeExact(chain,branch,dnawidth = 2,helang=135):
    # place beads and branches exactly as specified by the arrays (no interpolation)
    atoms1 = []; atoms2 = []; atoms3= []
    count = 0
    theta = helang*pi/180
    print len(chain)
    for bp in range(len(chain)):
        count += 1
        atoms3.append(Atom(coords=chain[bp][:],num=count,name='A3'))
        if len(atoms3) > 1:
            atoms3[-1].conect.append(atoms3[-2])
            atoms3[-2].conect.append(atoms3[-1])

            if bp < len(chain)-1:
                tvec = [chain[bp+1][i] - chain[bp][i] for i in range(3)]
            else:
                tvec = [chain[bp][i] - chain[bp-1][i] for i in range(3)]
            # x axis at this particular basepair
            xax = [branch[bp][i] - chain[bp][i] for i in range(3)]
            yax = crossprod(tvec,xax)
            yax = normalize(yax)
            coords = [chain[bp][i] + xax[i] for i in range(3)]
            count += 1
            atoms1.append(Atom(coords=coords,num=count,name='A1'))

            count += 1
            coords = [chain[bp][i] + dnawidth/2*(cos(theta)*xax[i] + sin(theta)*yax[i]) for i in range(3)]
            atoms2.append(Atom(coords=coords,num=count,name='A2'))

            atoms1[-1].resnum = bp+1
            atoms2[-1].resnum = bp+1
            atoms3[-1].resnum = bp+1         

            # keep track of atom connections
            atoms1[-1].conect.append(atoms2[-1])
            atoms2[-1].conect.append(atoms1[-1])
            if len(atoms1) > 1:
                atoms1[-1].conect.append(atoms1[-2])
                atoms1[-2].conect.append(atoms1[-1])
            if len(atoms2) > 1:
                atoms2[-1].conect.append(atoms2[-2])
                atoms2[-2].conect.append(atoms2[-1])

    return atoms1+atoms2+atoms3

def getNucAtoms(bends,nuch=5.45):
    # get the atoms for defining the nucleosomes
    # based on the center positions and the quaternion orientations

    centAtm = []; cylAtm1 = []; cylAtm2 =[]; tailAtm = []; chrgAtm = []; mptAtm = []

    for bnd in bends:    
        centAtm.append(Atom(coords=bnd.cent,num=1,name='AC'))

        # get z axis from quaternions
        a = bnd.orient[0]; b = bnd.orient[1]; c = bnd.orient[2]; d = bnd.orient[3]
        z = [2*a*c+2*b*d,2*c*d - 2*a*b , a*a-b*b-c*c+d*d]
        nz = norm(z)
        z = [z[i]/nz for i in range(3)]        
        atm1 = [bnd.cent[i]+nuch/2*z[i] for i in range(3)];
        atm2 = [bnd.cent[i]-nuch/2*z[i] for i in range(3)];
        cylAtm1.append(Atom(coords=atm1,num=1,name='AC1'))
        cylAtm2.append(Atom(coords=atm2,num=1,name='AC2'))
        
        # connections
        centAtm[-1].conect.append(cylAtm1[-1])
        centAtm[-1].conect.append(cylAtm2[-1])
        cylAtm1[-1].conect.append(centAtm[-1])
        cylAtm2[-1].conect.append(centAtm[-1])

        # get the flexible tail atoms
        for tail in bnd.tails:
            for c in range(len(tail)):
                tailAtm.append(Atom(coords=tail[c],num=1,name='FT'))
                if c == 0:
                    #centAtm[-1].conect.append(tailAtm[-1])
                    #tailAtm[-1].conect.append(centAtm[-1])
                    pass
                else:
                    tailAtm[-1].conect.append(tailAtm[-2])
                    tailAtm[-2].conect.append(tailAtm[-1])
                    
        # get the charge atoms
        for charge in bnd.charges:
            chrgAtm.append(Atom(coords=charge,num=1,name='CH'))
              
        # get the marked point atoms
        for pt in bnd.markpoints:
            mptAtm.append(Atom(coords=pt[1],num=1,name=pt[0]))

    return [centAtm, cylAtm1, cylAtm2, tailAtm, chrgAtm,mptAtm]

def getCoords(infile):
    # get the chain and branch coordinates from a file
    # also return nucleosome centers and quaternions
    IF = open(infile)

    chain = []; branch = []; bends = []; atmtypes = []; quaterns=[];
    while 1:
        line = IF.readline()
        if len(line) == 0:
            break
        data = line.split()
        if data[0] == 'A' or data[0] == 'F':
            chain.append([float(i) for i in data[1:4]])
            branch.append([float(i) for i in data[4:]])
            atmtypes.append(data[0])
        elif data[0] == 'B':
            # read in bend center and orientation
            bends.append(Bend([float(i) for i in data[1:4]],[float(i) for i in data[4:8]]))
        elif data[0] == 'T':
            # read in tail bend
            if int(data[-1])==0:
                # this is a new tail
                bends[-1].tails.append([( float(data[1]),float(data[2]),float(data[3]) )])
            else:
                bends[-1].tails[-1].append(( float(data[1]),float(data[2]),float(data[3]) ))
        elif data[0] == 'C':
            # read in coordinates for nucleosome charges 
            bends[-1].charges.append([float(data[1]),float(data[2]),float(data[3])])
        elif data[0][0] == 'M':
            # other marked points
            bends[-1].markpoints.append([data[0][1:],[float(data[1]),float(data[2]),float(data[3])]])
        else:
            print 'ERROR: type of line must be A,B,F,T,C, or start with M', data
            sys.exit()
            
    IF.close()
    return [chain,branch,atmtypes,bends]

# ---------------------------------------
# -------- run the actual code ----------
# ---------------------------------------

def runcode(argv):
	if len(argv) < 2:
	    print """usage: beads2pdb.py infile(s) -o <outfile> -nbp <bp> -fnbp <fnbp> -nucfile <nucfile> \
-scl <scl> -twist <twist> -nucfile <nucfile.pdb> -nr -connectmax <maxdist>
            All argument pairs besides the input file are optional. Can also supply a list or glob (eg: file.*.out) of input files. Input files may not start with "-"

            -o Gives the output pdb file name. Default output file is: infileroot.pdb

	    -nbp number of base pairs per linker DNA segment. Default is 1.

	    -fnbp number of base pairs per filler (bound) DNA segment. Default is 1

	    -twist Expected DNA twist density, which is used for figuring out in 
             which direction to interpolate branches. Default is 2*pi/(10.46*0.34)

            -nuch Height of the nucleosome axis (default 5.5)

            -scl is the scaling factor to convert from .out to .pdb coordinates 

             (default 10, as the .out file is expected to be in nm while the .pdb file is in angstroms)

            -nucfile Place actual nucleosome coordinates, as read in from a pdb file,
             at the position of the nucleosomes in the structure

            -nr turn off the recentering of the coordinates to shift the COM at the origin

            -connectmax is the maximal distance such that atoms that would normally be 
             connected are not if the distance betweent them is greater than this. 
             This is useful for single-nucleosome calculations to prevent connection 
             between the linker DNA ends and the distant edge filler beads"""            
	    sys.exit()

        infiles = []
        for a in sys.argv[1:]:
            if a[0]=='-':
                break
            infiles.append(a)           
        
	if '-o' in argv:
	    ind = argv.index('-o')
	    outfiles = [argv[ind+1]]*len(infiles)
	else:
            outfiles = infiles[:]
            for c in range(len(infiles)):
                (root,ext) = os.path.splitext(infiles[c])
                outfiles[c] = root+'.pdb'
        
	if '-nbp' in argv:
	    ind = argv.index('-nbp')
	    nbp = int(argv[ind+1])
	else:
	    nbp = 1

	if '-fnbp' in argv:
	    ind = argv.index('-fnbp')
	    fillnbp = int(argv[ind+1])
	else:
	    fillnbp = 1

	if '-twist' in argv:
	    ind = argv.index('-twist')
	    twist = float(argv[ind+1])
	else:
	    twist = 2*pi/(10.46*0.34)

        if '-nr' in argv:
            recenter = 0
        else:
            recenter = 1

        if '-scl' in argv:
            ind = argv.index('-scl')
            scl = argv[ind+1]
        else:
            scl = 10

        if '-nuch' in argv:
            ind = argv.index('-nuch')
            nuch = float(argv[ind+1])
        else:
            nuch = 5.5
        
        if '-connectmax' in argv:
            ind = argv.index('-connectmax')
            connectmax = float(argv[ind+1])
        else:
            connectmax = -1

        placenuc = '-nucfile' in argv
        if placenuc:
            ind = argv.index('-nucfile')
            nucpdbfile = argv[ind+1]

        for cf in range(len(infiles)):
            infile = infiles[cf]; outfile = outfiles[cf]
            print "Input file, output file:", infile, outfile

            [chain,branch,atmtypes,bends] = getCoords(infile)              

	    # place all the DNA atoms
            if '-exact' in argv:
                atoms = placeExact(chain,branch)
            else:
                atoms = placeBP(chain,branch,atmtypes,nbp=nbp,fillnbp=fillnbp,twist=twist,connectmax=connectmax)                      

            [centatms,cylAtm1,cylAtm2,tailAtm,chrgAtm,mptAtm] = getNucAtoms(bends,nuch)	
            atoms.extend(centatms); atoms.extend(cylAtm1)
            atoms.extend(cylAtm2); atoms.extend(tailAtm)
            atoms.extend(chrgAtm)
            atoms.extend(mptAtm)

            if '-orientnucz' in argv:
                # reorient everything so that 1st nucleosome axis is along z axis
                # and initial bead is along x-axis 
                cent = centatms[0].coords  
 
                zax = cylAtm2[0].coords-cent
                zax = zax/sqrt(dot(zax,zax))
            
                centatms = [a for a in atoms if a.name=='A3']
                xax = centatms[0].coords - cent
                xax = xax - dot(xax,zax)*zax
                xax = xax/sqrt(dot(xax,xax))
                yax = cross(zax,xax)

                # rotation matrix to invert to canonical coords
                rotmat = array([xax,yax,zax])

                for a in atoms:
                    tmp = a.coords - cent
                    a.coords = dot(rotmat,tmp)                

            if (scl != 1):
                for a in atoms:
                    a.coords = a.coords*scl

            for a in atoms:
                a.type = 'ATOM'

            fullstruct = makeBareStruct()
            #fullstruct.atoms = atoms               

            # place the actual crystal structure nucleosome coordinates
            if placenuc:
                nucstruct = Structure(nucpdbfile)                    

                # get rid of extra atoms in nucleosome struct
                nucstruct.atoms = [a for a in nucstruct.atoms if not( a.type == 'HETATM' and a.name not in ['PC','MC'])]
                nucstruct.recenter()

                for b in range(len(centatms)):
                    # rotation matrix from quaternion
                    q= Quaternion(bends[b].orient)
                    mat = q.getRotMat()
                    mat = matrix(mat).T
                    for a in nucstruct.atoms:                   
                        atm = copy(a)
                        # shift and rotate as appropriate
                        atm.coords = array(atm.coords*mat)[0]
                        atm.coords = atm.coords+centatms[b].coords
                        fullstruct.atoms.append(atm)

                fullstruct.startlines = nucstruct.startlines

            # renumber the atoms
            fullstruct.atoms.extend(atoms)
            fullstruct.renumAtoms()

	    # recenter the structure
            if recenter:
                fullstruct.recenter()
              
            # output the PDB file
                fullstruct.outputPDB(outfile)

if __name__ == '__main__':
	runcode(sys.argv)

