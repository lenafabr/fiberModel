""" generic utilities for python scripts """
from numpy import *
from math import *

class Quaternion:
    def __init__(self,coords):
        self.coords = array(coords)
    
    def getRotMat(self):
        """get rotation matrix for this quaternion"""
        w = self.coords[0]
        x = self.coords[1]
        y = self.coords[2]
        z = self.coords[3];

        mat = [[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(w*y + x*z)], \
                   [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)], \
                   [2*(x*z-w*y), 2*(w*x+y*z), 1 - 2*(x*x+y*y)]]

        return array(mat)

def norm(x):
    """ get the norm of x """
    if type(x).__name__ == 'ndarray':
        return sqrt(dot(x,x))
    else:
        ans = 0
        for i in x:
            ans += i**2
        ans = sqrt(ans)
        return ans

def normalize(x):
    """get norm of x and divide through by it"""
    n = norm(x)
    x = [i/n for i in x]
    return x

def dotprod(x,y):
    """dot product of two vectors"""

    if type(x).__name__ == 'ndarray':
        return dot(x,y)
    else:
        return sum([x[i]*y[i] for i in range(len(x))])

def crossprod(a,b):
    """returns cross product of 2 triplets"""
    
    x = a[1]*b[2] - a[2]*b[1]
    y = a[2]*b[0] - a[0]*b[2]
    z = a[0]*b[1] - a[1]*b[0]

    return [x,y,z]

def dihedral(i,j,k,l):
    """dihedral angle between coordinates i-j-k-k"""

    b1 = [j[a]-i[a] for a in range(3)]
    b2 = [k[a]-j[a] for a in range(3)]
    b3 = [l[a]-k[a] for a in range(3)]

    phi = atan2(norm(b2)*dotprod(b1,crossprod(b2,b3)),dotprod(crossprod(b1,b2),crossprod(b2,b3)))

    return phi

def mat2euler(mat):
    """Convert from a rotation matrix to z-x-z euler angles"""

    eul = array([0.0,0.0,0.0])
    if abs(mat[2,2])==1:
        eul[0] = atan2(mat[1,0],mat[0,0])
    else:
        eul[0] = atan2(mat[0,2],-mat[1,2])
        eul[1] = acos(mat[2,2])
        eul[2] = atan2(mat[2,0],mat[2,1])

    return eul

def rotmat2AngleAxis(mat):
    # convert from a rotation matrix to an angle-axis representation
    # from Diebel, except my matrix is the transpose of theirs
    r11 = mat[0,0]; r12 = mat[1,0]; r13 = mat[2,0]
    r21 = mat[0,1]; r22=mat[1,1]; r23 = mat[2,1];
    r31 = mat[0,2]; r32 = mat[1,2]; r33 = mat[2,2]

    if r22 > -r33 and r11 > -r22 and r11>-r33:
        tmp = sqrt(1.0+r11+r22+r33)
        qvec = 0.5*array([tmp,(r23-r32)/tmp, (r31-r13)/tmp,(r12-r21)/tmp])     
    elif r22<-r33 and r11>r22 and r11>r33:
        tmp = sqrt(1.0+r11-r22-r33)
        qvec = 0.5*array([(r23-r32)/tmp,tmp,(r12+r21)/tmp,(r31+r13)/tmp])
    elif r22>r33 and r11<r22 and r11<-r33:
        tmp = sqrt(1.0-r11+r22-r33)
        qvec = 0.5*array([(r31-r13)/tmp,(r12+r21)/tmp,tmp,(r23+r32)/tmp])
    elif r22<r33 and r11<-r22 and r11<r33:
        tmp = sqrt(1.0-r11-r22+r33)
        qvec = 0.5*array([(r12-r21)/tmp,(r31+r13)/tmp,(r23+r32)/tmp,tmp])
    else:
        print "Bad rotation matrix:\n%s" %mat
                    
    angle = acos(qvec[0])*2
    axis = qvec[1:4]/norm(qvec[1:4])

    return [angle,axis]
