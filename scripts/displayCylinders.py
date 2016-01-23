# display the nucleosomes as cylinders
import re
from pymol.cgo import *
from pymol import cmd, stored
from math import *

atomR = re.compile('^(ATOM|HETATM) *([0-9]+) *(\S+) +([A-Z]+) +([A-Z]?) +\
([-0-9]+) *(-?[0-9.]+) *(-?[0-9.]+) *(-?[0-9.]+)(.*$)')

def getCylinderCoords(cylinderCenterName, cylinderEdgeNames):
	"""Queries pymol for the specified coordinates representing the cylinders.
	
	cylinderCenterName is a string specifying the pdb atom name field of the cylinder centers
	cylinderEdgeNames is a two element list specifying the pdb atom names of the cylinder top and bottom
	
	returns nothing, stores the results in pymol.stored.centCoords, pymol.stored.cylCoords0, pymol.stored.cylCoords1
	"""
	
	
	stored.centCoords = [];
	stored.cylCoords0 = [];
	stored.cylCoords1 = [];
	
	cyl0Sel = "cylAtoms0";
	cyl1Sel = "cylAtoms1";
	centerSel = "centers";
	
	cmd.select(centerSel, "(name " + cylinderCenterName + ")");
	cmd.select(cyl0Sel, "(name " + cylinderEdgeNames[0] + ")");
	cmd.select(cyl1Sel, "(name " + cylinderEdgeNames[1] + ")");
	
	cmd.iterate_state(1, selector.process(centerSel), "stored.centCoords.append([x,y,z])");
	cmd.iterate_state(1, selector.process(cyl0Sel), "stored.cylCoords0.append([x,y,z])");
	cmd.iterate_state(1, selector.process(cyl1Sel), "stored.cylCoords1.append([x,y,z])");


def makeCylinders(cylinderCenterName, cylinderEdgeName0, cylinderEdgeName1, rad=41.8, height=None, transparency=0):
	"""Makes cylindrical representations of structures based on the specified coordinates.
	
	cylinderCenterName is a string specifying the pdb atom name field of the cylinder centers
	cylinderEdgeName0 and cylinderEdgeName1 are strings specifying the pdb atom name field of the cylinder top and bottom
	rad is the radius of the cylinder in angstroms (optional, default: 41.8)
	height is the height of the cylinder in angstroms (optional, default: determined by cylinder edge coords)

	adds the cylinders to the pymol display with the name 'cylinders'

	"""

	rad = float(rad); transparency = float(transparency);

	getCylinderCoords(cylinderCenterName, [cylinderEdgeName0, cylinderEdgeName1]);
	
	obj = [];
	
	for i in range(len(stored.centCoords)):
		
		cent = stored.centCoords[i];
		cyl0 = stored.cylCoords0[i];
		cyl1 = stored.cylCoords1[i];
		
		if (not (height == None)):
			height = float(height);
			zaxis = [cyl0[j] - cent[j] for j in range(len(cyl0))];
			nz = sqrt(zaxis[0]**2+zaxis[1]**2+zaxis[2]**2);
			zaxis = [i/nz for i in zaxis];
			cyl0 = [cent[j] + zaxis[j] * height/2.0 for j in range(len(cent))];
			cyl1 = [cent[j] - zaxis[j] * height/2.0 for j in range(len(cent))];
			
		obj.extend([ALPHA,1-transparency]);
		nextCyl = [CYLINDER, cyl0[0], cyl0[1], cyl0[2], cyl1[0], cyl1[1], cyl1[2], rad, 0.96, 0.87, 0.70, 0.96,0.87,0.70];
		obj.extend(nextCyl);
		
		
	cmd.load_cgo(obj, 'cylinders');

cmd.extend("makeCylinders", makeCylinders);
