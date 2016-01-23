#!/usr/bin/python
# Run test calculations and check results
import sys
try:
    import subprocess
except ImportError:
    print "Your version of python does not include the subprocess module. You will need a newer Python version (2.6.* recommended)"
    sys.exit()

tiny = 1e-4
small = 1e-2

class TestError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def readOptimFile(infile):
    """Read in an output file from an OPTIMIZE or GETSTRUCT calculation. Returns dictionary with energy, coordinates, diameter, energy parts, nuc dists, initial energy if present"""
    
    IF = open(infile, 'r')
    lines = IF.readlines()
    IF.close()

    datadict = {}

    prevline = ''
    for line in lines:
        if line.strip().startswith('Chain energy'):
            data = line.split()
            try:
                datadict['energy'] = float(data[-1])
            except:
                print "Problem reading chain energy in output file: "
                print line                    
                raise IOError()
            
        if line.strip().startswith('Overall fiber diameter'):
            data = line.split()
            try:
                datadict['diam'] = float(data[-1])
            except:
                print "Problem reading fiber diameter in output file: "
                print line                    
                raise IOError()
           
        elif line.strip().startswith('Regular helix coordinates'):
            datadict['helcrds'] = [float(i) for i in line.split()[-6:]]

        elif line.strip().startswith('Energy parts'):
            datadict['eparts'] = [float(i) for i in line.split()[-6:]]

        elif prevline.strip().startswith('Nucleosome distances'):
            datadict['nucdist'] = [float(i) for i in line.split()[-20:]]

        elif line.strip().startswith('Initial chain energy'):
            datadict['initenergy'] = float(line.split()[-1])

        elif line.strip().startswith('Final energy'):
            datadict['energy'] = float(line.split()[-1])

        elif line.strip().startswith('Final helix parameters'):
            datadict['helcrds'] = [float(x) for x in line.split()[-6:]]

        prevline = line

    return datadict

def readDatabaseFile(infile):
    """Read in an output file from a database parsing or basinhopping calculation. Return results as a dictionary of values, for all final configurations"""

    IF = open(infile,'r')
    lines = IF.readlines()
    IF.close()
    
    datadict={}
    firsthop = 1
    for c in range(len(lines)):
        line = lines[c]

        # starting hop
        if line.strip().startswith('HOP, MINE') and firsthop:
            datadict['starthop'] = int(line.split()[-7])           
            firsthop = 0
        
        # energies of structures
        if line.strip().startswith('CHAIN,'):
            j = c+2
            datadict['energy'] = []
            while not '-----------' in lines[j]:
                datadict['energy'].append(float(lines[j].split()[-1]))
                j+=1
            datadict['nstruct'] = len(datadict['energy'])            

        # helix coordinates 
        if line.strip().startswith('HEIGHT,'):
            j = c+2
            datadict['helcrds'] = []
            while not '-----------' in lines[j]:
                datadict['helcrds'].append([float(x) for x in lines[j].split()[-7:]])
                j+=1

        # energy parts
        if line.strip().startswith('BENDING,'):
            j = c+2
            datadict['eparts'] = []
            while not '-----------' in lines[j]:
                datadict['eparts'].append([float(x) for x in lines[j].split()[-6:]])
                j+=1

        # nucleosome distances
        if line.strip().startswith('NUCLEOSOME DISTANCES'):
            j = c+2
            datadict['nucdist'] = []
            while not '-----------' in lines[j]:
                datadict['nucdist'].append([float(x) for x in lines[j].split()[1:]])
                j+=1       

    return datadict

def runCalculation(suf,exe='../fibermodel.exe'):
    """Run a particular calculation with the given suffix. Check return code and stderr file"""
    
    OF = open('stdout.%s' %suf,'w')
    EF = open('stderr.%s' %suf,'w')
    proc = subprocess.Popen([exe,suf], stdout=OF,stderr=EF)
    proc.wait()        
    OF.close()
    EF.close()
    
    if proc.returncode != 0:
        val =  "Bad return code: %d" %proc.returncode
        return TestError(val)  

    # check error file
    EF = open('stderr.%s' %suf, 'r')    
    errtext = EF.read()        
    EF.close()
    if len(errtext)>0:
        val =  "Some errors found in error file: \n%s" %errtext
        raise TestError(val)

    return 0

def straightTest(suf):
    """Test straight-linker structure"""      

    print "Running test of straight-linker structure (this may take a couple minutes)..."

    # run calculation
    runCalculation(suf)

    # compare output to the reference file
    refdict = readOptimFile('stdout.%s.ref' %suf)
    datadict = readOptimFile('stdout.%s' %suf)

    # compare energies
    if abs(refdict['energy']-datadict['energy']) > tiny:
        raise TestError("Wrong energy calculated. Found %20.10f. Should be %20.10f" %(datadict['energy'], refdict['energy']))

    # compare energy parts
    for c in range(6):
        if abs(refdict['eparts'][c] - datadict['eparts'][c]) > tiny:
            raise TestError("Wrong energy part %d. Found %20.10f. Should be %20.10f" %(c+1, datadict['eparts'][c], refdict['eparts'][c]))

    # compare helix coordinates
    for c in range(6):
        if abs(refdict['helcrds'][c] - datadict['helcrds'][c]) > tiny:
            raise TestError("Wrong helix coord %d. Found %20.10f. Should be %20.10f" %(c+1, datadict['helcrds'][c], refdict['helcrds'][c]))

    # compare fiber diameter
    if abs(refdict['diam']-datadict['diam']) > tiny:
        raise TestError("Wrong fiber diameter calculated. Found %20.10f. Should be %20.10f" %(datadict['diam'], refdict['diam']))

    return 0
                     
def energyTest(suf):
    """Test energy calculation for a restarted structure"""
    
    print "Running test of energy calculation..."

    # run calculation
    runCalculation(suf)

    # compare output to the reference file
    refdict = readOptimFile('stdout.%s.ref' %suf)
    datadict = readOptimFile('stdout.%s' %suf)

    # compare energies
    if abs(refdict['energy']-datadict['energy']) > tiny:
        raise TestError("Wrong energy calculated. Found %20.10f. Should be %20.10f" %(datadict['energy'], refdict['energy']))

    # compare energy parts
    for c in range(6):
        if abs(refdict['eparts'][c] - datadict['eparts'][c]) > tiny:
            raise TestError("Wrong energy part %d. Found %20.10f. Should be %20.10f" %(c+1, datadict['eparts'][c], refdict['eparts'][c]))

    # compare helix coordinates
    for c in range(6):
        if abs(refdict['helcrds'][c] - datadict['helcrds'][c]) > tiny:
            raise TestError("Wrong helix coord %d. Found %20.10f. Should be %20.10f" %(c+1, datadict['helcrds'][c], refdict['helcrds'][c]))

    # compare fiber diameter
    if abs(refdict['diam']-datadict['diam']) > tiny:
        raise TestError("Wrong fiber diameter calculated. Found %20.10f. Should be %20.10f" %(datadict['diam'], refdict['diam']))

    # compare nucleosome distances
    for c in range(len(refdict['nucdist'])):
        if abs(refdict['nucdist'][c] - datadict['nucdist'][c]) > tiny:
            raise TestError("Wrong nucleosome distance %d. Found %20.10f. Should be %20.10f" %(c+1, datadict['nucdist'][c], refdict['nucdist'][c]))

    return 0
    

def optimTest(suf):
    """Test optimization calculation"""    

    print "Running test of optimization..."

    # run calculation
    runCalculation(suf)

    # compare output to the reference file
    refdict = readOptimFile('stdout.%s.ref' %suf)
    datadict = readOptimFile('stdout.%s' %suf)

    # compare initial energies
    if abs(refdict['initenergy']-datadict['initenergy']) > tiny:
        raise TestError("Wrong initial energy. Found %20.10f. Should be %20.10f" %(datadict['initenergy'], refdict['initenergy']))

    # compare energies
    if abs(refdict['energy']-datadict['energy']) > small:
        raise TestError("Wrong energy calculated. Found %20.10f. Should be %20.10f" %(datadict['energy'], refdict['energy']))

    # compare energy parts
    for c in range(6):
        if abs(refdict['eparts'][c] - datadict['eparts'][c]) > small:
            raise TestError("Wrong energy part %d. Found %20.10f. Should be %20.10f" %(c+1, datadict['eparts'][c], refdict['eparts'][c]))

    # compare helix coordinates
    for c in range(6):
        if abs(refdict['helcrds'][c] - datadict['helcrds'][c]) > small:
            raise TestError("Wrong helix coord %d. Found %20.10f. Should be %20.10f" %(c+1, datadict['helcrds'][c], refdict['helcrds'][c]))

    # compare fiber diameter
    if abs(refdict['diam']-datadict['diam']) > small:
        raise TestError("Wrong fiber diameter calculated. Found %20.10f. Should be %20.10f" %(datadict['diam'], refdict['diam']))   

    return 0

def basinhopTest(suf):
    """Test basin-hopping calculation"""
    
    print "Running basin-hopping test..."
    
    # run calculation
    runCalculation(suf) 

    # compare results
    refdict = readDatabaseFile('stdout.%s.ref' %suf)
    datadict = readDatabaseFile('stdout.%s' %suf)

    # check initial hop
    if refdict['starthop'] != datadict['starthop']:
        raise TestError("Started with the wrong hop %d. Should be %d" %(datadict['starthop'], refdict['starthop']))

    # which new structures correspond to the first 2 reference structures?
    structind = [-1,-1]
    for i in range(2):
        for c in range(datadict['nstruct']):
            if abs(datadict['energy'][c]-refdict['energy'][i]) < small:
                structind[i] = c
                break
        if structind[i] < 0:
            raise TestError('Failed to find a structure with energy %20.10f corresponding to reference structure %d' %(refdict['energy'][i], i+1))

    print "Corresponding structures:", structind

    for i in range(2):
        j = structind[i]

        # check energy parts
        for c in range(6):
            if abs(refdict['eparts'][i][c] - datadict['eparts'][j][c]) > small:
                raise TestError("Wrong energy part %d in structure %d. Found %20.10f. Should be %20.10f" %(c+1, i+1, datadict['eparts'][j][c], refdict['eparts'][i][c]))

        # check helix coords
        for c in range(3)+[4,6]:
            if abs(refdict['helcrds'][i][c] - datadict['helcrds'][j][c]) > small:
                raise TestError("Wrong helix coord %d in structure %d. Found %20.10f. Should be %20.10f" %(c+1, i+1, datadict['helcrds'][j][c], refdict['helcrds'][i][c]))

    return 0

def databaseTest(suf):
    """Test database parsing calculation"""
    
    print "Running database parsing test..."
    
    # run calculation
    runCalculation(suf) 

    # compare results
    refdict = readDatabaseFile('stdout.%s.ref' %suf)
    datadict = readDatabaseFile('stdout.%s' %suf)
   
    # compare first 3 structures
    for i in range(3):
        # check energies
        if abs(refdict['energy'][i] - datadict['energy'][i]) > small:
            raise TestError("Wrong energy for structure %d. Found %20.10f. Should be %20.10f" %(i+1,datadict['energy'][i], refdict['energy'][i]))

        # check energy parts
        for c in range(6):
            if abs(refdict['eparts'][i][c] - datadict['eparts'][i][c]) > small:
                raise TestError("Wrong energy part %d in structure %d. Found %20.10f. Should be %20.10f" %(c+1, i+1, datadict['eparts'][i][c], refdict['eparts'][i][c]))

        # check helix coords
        for c in range(3)+[4,6]:
            if abs(refdict['helcrds'][i][c] - datadict['helcrds'][i][c]) > small:
                raise TestError("Wrong helix coord %d in structure %d. Found %20.10f. Should be %20.10f" %(c+1, i+1, datadict['helcrds'][i][c], refdict['helcrds'][i][c]))

        # check nucleosome distances
        for c in range(len(refdict['nucdist'][i])):
            if abs(refdict['nucdist'][i][c] - datadict['nucdist'][i][c]) > small:
                raise TestError("Wrong nucleosome distance %d in structure %d. Found %20.10f. Should be %20.10f" %(c+1, i+1, datadict['nucdist'][i][c], refdict['nucdist'][i][c]))

    return 0

# -----------------------
# run the tests 
# ----------------------

allgood = 1
                        
try:
    res = straightTest('teststraight')
    print "Straight-linker test successful!"
except TestError, e:
    allgood = 0
    print "Failed straight-linker test! Error message: %s" %e.value
    print "See stdout.teststraight and stderr.teststraight for more info"

try:
    res = energyTest('testenergy')
    print "Energy calculation test successful!"
except TestError, e:
    allgood = 0
    print "Failed energy calculation test! Error message: %s" %e.value
    print "See stdout.testenergy and stderr.testenergy for more info"

try:
    res = optimTest('testoptim')
    print "Optimization test successful!"
except TestError, e:
    allgood = 0
    print "Failed optimization test! Error message: %s" %e.value
    print "See stdout.testoptim and stderr.testoptim for more info"

try:
    res = basinhopTest('testbasinhop')
    print "Basin hopping test successful!"
except TestError, e:
    allgood = 0
    print "Failed basin hopping test! Error message: %s" %e.value
    print "See stdout.testbasinhop and stderr.testbasinhop for more info"

try:
    res = databaseTest('testdbparse')
    print "Database parsing test successful!"
except TestError, e:
    allgood = 0
    print "Failed database parsing test! Error message: %s" %e.value
    print "See stdout.testdbparse and stderr.testdbparse for more info"

if allgood:
    print "All tests successfully completed. The fiberModel code works properly on this system."
else:
    print "One of the tests failed. Something is not quite right with the code, or it might just be a slight numerical error. Use at own risk."
