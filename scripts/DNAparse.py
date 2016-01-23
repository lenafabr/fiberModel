# functions for extracting coarse-grained parameters (such as axis an orientation) from DNA coordinates in a pdb file

from pdbutils import *
import copy
from genutils import *

class DNAobj(Structure):
    """Object corresponding to a DNA double-helix"""
    def __init__(self):        
        self.initvars()
        self.bpairs = []
        self.pairind = []

    def DNAFromStruct(self,struct,DNAres = ['DA','DT','DC','DG']):
        """ Split out just the DNA chains from a structure object. DNAres specifies which residues count as DNA"""
        
        self.chains = [];  self.atoms = []
        self.residues = [r for r in struct.residues if r.name in DNAres]
        self.resetFromResidues()

        # self.chains = []; self.residues = []; self.atoms = []
        # for c in struct.chains:
        #     if c.residues[0].name in DNAres:
        #         self.chains.append(c)

        # for c in self.chains:
        #     self.residues.extend(c.residues)
        #     for r in self.residues:
        #         self.atoms.extend(r.atoms)

        self.startlines = struct.startlines
        self.endlines = struct.endlines

        return 0
        
    def setup(self,infile,DNAonly = 0):
        # given an input file, set up the DNA molecule based on the chains found in that file
        # for now, will fail  only use first 2 chains found
        # if DNAonly=1, then assume all atoms in the file are from DNA (so residues don't need to start with D)

        self.structFromFile(infile)

        if (len(self.chains) < 2):
            raise IOError("Error! Could not find 2 DNA chains!")
        
        if not DNAonly: # filter out all the non DNA stuff            
            self.residues = []; self.atoms = []
            # pick out only the DNA residues
            for c in self.chains:
                c.residues = [r for r in c.residues if r.name[0] == 'D']
                self.residues.extend(c.residues)
                for r in c.residues:
                    self.atoms.extend(r.atoms)

            # renumber the atoms
            for i in range(len(self.atoms)):
                self.atoms[i].num = i+1
            
        self.chains = [c for c in self.chains if len(c.residues) > 0]
        
        return 0

    def setPairs(self,reverse=None):
        # set the basepairs
        # if reverse=true, then assume chains are both given 3' to 5' (or vice versa)
        # if reverse is None, then decide whether to reverse or not based on sequence

        nbp = len(self.chains[0].residues)
        if reverse == None:
            dorev = 0
            for a in range(nbp):
                if not goodPair(self.chains[0].residues[a], self.chains[1].residues[a]):
                    dorev = 1
                    break
        else:
            dorev = reverse

        if dorev:
            for a in range(nbp):
                self.pairind.append(nbp-a-1)
        else:
            for a in range(nbp):
                self.pairind.append(a)

        for i in range(nbp):
            r1 = self.chains[0].residues[i]; r2 = self.chains[1].residues[self.pairind[i]]
            if not goodPair(r1,r2):                
                print "WARNING: base pairs don't match up: ", r1.name, r2.name

            bp = BasePair([self.chains[0].residues[i],self.chains[1].residues[self.pairind[i]]])
            self.bpairs.append(bp)

        return 0 

    def getInmostN(self):
        # for all the basepairs pick out the inmost N atoms
        for b in self.bpairs:
            b.inN = []
            for r in b.residues:
                if r.name[-1] in ['A' ,'G']:
                    b.inN.append(r.getAtomByName('N1'))
                else:
                    b.inN.append(r.getAtomByName('N3'))
                if b.inN[-1] == -1:
                    raise ValueError('Failed to find inmost nitrogen for residues %s' %r)

        return 0

    def getAxis(self):
        # get the points defining the center axis of the DNA, one point for each basepair
        if len(self.bpairs) == 0:            
            self.setPairs()
        if not hasattr(self.bpairs[0],'inN'):
            self.getInmostN()

        self.axis = [(b.inN[0].coords+b.inN[1].coords)/2 for b in self.bpairs]

        # get atoms corresponding to the DNA axis
        axatoms = []; num = self.atoms[-1].num
        for i in range(len(self.axis)):
            newatom = Atom(coords=self.axis[i])
            num = num + 1
            newatom.num = num
            axatoms.append(newatom)
            if i > 0:
                newatom.conect.append(axatoms[i-1])                
                axatoms[i-1].conect.append(axatoms[i])       
            axatoms[i].name = 'AX'
        
        return axatoms
  
    def getTangents(self,offset=0,bpback=1,axlen=10):
        """Get the tangents and the x-axes at the DNA ends. Offset is how many base pairs to ignore at the edges. bpback is how many basepairs to reach back in calculating the tangent vector. Minus end tangent points inward towards the bend.axlen gives the length of the axes for placing the checking atoms"""

        tan1 = self.axis[offset+bpback] - self.axis[offset]
        tan1 = tan1/norm(tan1)
        print 'Minus end tangent: ', tan1

        centatm1 = Atom(coords=self.axis[offset])        
        tanatm1 = Atom(coords=self.axis[offset]+tan1*axlen)
        tanatoms = [centatm1, tanatm1]
        centatm1.name = 'XC'
        tanatm1.name = 'XT'
        centatm1.conect = [tanatm1]
        tanatm1.conect =[centatm1]

        tan2 = self.axis[-1-offset] - self.axis[-1-offset-bpback]
        tan2 = tan2/norm(tan2)
        print 'Plus end tangent: ', tan2

        centatm2 = Atom(coords=self.axis[-1-offset])
        tanatm2 = Atom(coords=self.axis[-1-offset]+tan2*axlen)
        tanatoms.extend([centatm2, tanatm2])
        centatm2.name = 'XC'
        tanatm2.name = 'XT'
        centatm2.conect = [tanatm2]
        tanatm2.conect =[centatm2]
      

        xatm1 = self.bpairs[offset].residues[0].getAtomByName("^C1'",1)
        xax1 = xatm1.coords - self.axis[offset]
        xax1 = xax1 - dot(xax1,tan1)*tan1
        xax1 = xax1/norm(xax1)
        print 'Minus end branch: ', xax1

        xatm1 = Atom(coords=centatm1.coords + axlen*xax1)
        tanatoms.append(xatm1)
        xatm1.name = 'XB'
        xatm1.conect.append(centatm1)
        
        xatm2 = self.bpairs[-1-offset].residues[0].getAtomByName("^C1'",1)
        xax2 = xatm2.coords - self.axis[-1-offset]
        xax2 = xax2 - dot(xax2,tan2)*tan2
        xax2 = xax2/norm(xax2)
        print 'Plus end branch: ', xax2

        xatm2 = Atom(coords=centatm2.coords+axlen*xax2)
        tanatoms.append(xatm2)
        xatm2.name = 'XB'
        xatm2.conect.append(centatm2)
        
        return [tanatoms, [tan1,xax1,tan2,xax2]]

    def getEndPos(self,com,offset=0):
        """Get positions of ends from the center of mass.Offset is the number of basepairs at each end to ignore."""
        comatm = Atom(coords=com,name='COM')        

        posm = self.axis[offset] - com
        posp = self.axis[-1-offset] - com
        #print 'POSM:', posm
        #print 'POSP:', posp

        return [comatm, posm,posp]

    def getLength(self,offset=0):
        """ Get the contour length of the DNA, excluding some number of basepairs from each end"""

        contlen = 0
        for i in range(offset+1,len(self.axis)-offset):
            contlen = contlen + norm(self.axis[i]-self.axis[i-1])

        self.length = contlen
        return contlen
    
# ------------- for testing only --------------
if __name__ == "__main__":
    print "Running test code:"
    struct = Structure()
    struct.structFromFile('./1KX5-nucleosome.pdb')
    #struct.structFromFile('./2bnaB.pdb')
    DNA = DNAobj()
    DNA.DNAFromStruct(struct)
    #DNA.DNAFromStruct(struct,['A','T','G','C'])

    com = struct.getCOM()
    axatoms = DNA.getAxis()
    [tanatoms,tanvecs] = DNA.getTangents(5)
    [comatm,posm,posp] = DNA.getEndPos(com)
    DNA.atoms.extend(axatoms)
    DNA.atoms.extend(tanatoms)
    DNA.atoms.append(comatm)
    DNA.renumAtoms()
    DNA.outputPDB('test.pdb')
