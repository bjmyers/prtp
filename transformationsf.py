import numpy as np
from numba import jit
import prtp.transformMod as tr
#'UniTuple(float64[:], 2)(float64[:],float64[:],float64[:],int32,int32)'
#@jit(cache=True,nopython=True)
def rotatevector(x,y,z,theta,axis):
    #This function rotates a vector in a right handed fashion
    #axis is 1,2,3 for x,y,z axis rotation
    
    if (axis==1):
        tempy = np.cos(theta)*y - np.sin(theta)*z
        z = np.sin(theta)*y + np.cos(theta)*z
        return x,tempy,z
    elif (axis==2):
        tempx = np.cos(theta)*x + np.sin(theta)*z
        z = -np.sin(theta)*x + np.cos(theta)*z
        return tempx,y,z
    else:
        tempx = np.cos(theta)*x - np.sin(theta)*y
        y = np.sin(theta)*x + np.cos(theta)*y
        return tempx,y,z


def rotateaxis(x,y,z,ux,uy,uz,theta):
    #This function rotates a vector in a right handed fashion
    #axis is given by input ux,uy,uz

    mag = np.sqrt(ux**2 + uy**2 + uz**2)
    with np.errstate(invalid='ignore'):
        ux /= mag
        uy /= mag
        uz /= mag
    
    c = np.cos(theta)
    s = np.sin(theta)
    tempx = (c+ux**2*(1-c))*x + (ux*uy*(1-c)-uz*s)*y + (ux*uz*(1-c)+uy*s)*z
    tempy = (uy*ux*(1-c)+uz*s)*x + (c+uy**2*(1-c))*y + (uy*uz*(1-c)-ux*s)*z
    z = (uz*ux*(1-c)-uy*s)*x + (uz*uy*(1-c)+ux*s)*y + (c+uz**2*(1-c))*z
    return tempx,tempy,z,ux,uy,uz

#'float[:,:](float64[:,:])'
#@jit(cache=True,nopython=True)  
def reflect(l,m,n,ux,uy,uz,ind=None):

    dot = ux*l + uy*m + uz*n
    l -= 2*dot*ux
    m -= 2*dot*uy
    n -= 2*dot*uz
    
    return l,m,n
    
#'UniTuple(float64[:], 2)(float64[:],float64[:],float64[:],int32,int32)'
def refract(l,m,n,ux,uy,uz,n1,n2):
    
    dot = ux*l + uy*m + uz*n
    negativedot = (dot<0)
    ux = np.where(negativedot,-ux,ux)
    uy = np.where(negativedot,-uy,uy)
    uz = np.where(negativedot,-uz,uz)
    dot = np.abs(dot)
    
    unfinishedrays = (dot != 1.0)
    t1 = np.where(unfinishedrays,np.arccos(dot),0)
    t2 = np.where(unfinishedrays,np.arcsin((n1/n2)*np.sin(t1)),0)
    cx = uy*n-m*uz
    cy = l*uz-ux*n
    cz = ux*m-l*uy
    
    l,m,n,cx,cy,cz = rotateaxis(l, m, n, cx, cy, cz, t2-t1)
    alpha = np.sqrt(l**2 + m**2 + n**2)
    l = np.where(unfinishedrays,l/alpha,l)
    m = np.where(unfinishedrays,m/alpha,m)
    n = np.where(unfinishedrays,n/alpha,n)
    
    return l,m,n,ux,uy,uz
    
#@jit(cache=True,nopython=True)
def transform(x,y,z,l,m,n,ux,uy,uz,tx,ty,tz,rx,ry,rz,coords=None):
    # Coordinate system transform, translations are done first, then rotations in XYZ order
    # Rotations act to rotate a surface via the right hand rule
    # tx,ty,tz are linear translations
    # rx,ry,rz are angular translations   

    x += tx
    y += ty
    z += tz
    
    x,y,z = rotatevector(x,y,z,rx,1)
    l,m,n = rotatevector(l,m,n,rx,1)
    ux,uy,uz = rotatevector(ux,uy,uz,rx,1)
    
    x,y,z = rotatevector(x,y,z,ry,2)
    l,m,n = rotatevector(l,m,n,ry,2)
    ux,uy,uz = rotatevector(ux,uy,uz,ry,2)
    
    x,y,z = rotatevector(x,y,z,rz,3)
    l,m,n = rotatevector(l,m,n,rz,3)
    ux,uy,uz = rotatevector(ux,uy,uz,rz,3)
    
    if coords is not None:
        #Define rotation and translation matrices
        rotm = rotationM(-rx,-ry,-rz)
        tranm = translationM(-tx,-ty,-tz)
        rotmi = rotationM(-rx,-ry,-rz,inverse=True)
        tranmi = translationM(tx,ty,tz)
        #Dot rotation into forward transform
        coords[0] = np.dot(rotm,coords[0])
        coords[1] = np.dot(np.dot(rotm,tranm),coords[1])
        coords[2] = np.dot(coords[2],rotmi)
        coords[3] = np.dot(coords[3],np.dot(tranmi,rotmi))
    
    return x, y, z, l, m, n, ux, uy, uz
    
    
def itransform(x,y,z,l,m,n,ux,uy,uz,tx,ty,tz,rx,ry,rz,coords=None):
    # This function inverts the transform function with the same arguments
    # transform(rays,0,10,0,2,0,0) will be undone by itransform(rays,0,10,0,2,0,0)
    
    x,y,z = rotatevector(x,y,z,rz,3)
    l,m,n = rotatevector(l,m,n,rz,3)
    ux,uy,uz = rotatevector(ux,uy,uz,rz,3)
    
    x,y,z = rotatevector(x,y,z,ry,2)
    l,m,n = rotatevector(l,m,n,ry,2)
    ux,uy,uz = rotatevector(ux,uy,uz,ry,2)
    
    x,y,z = rotatevector(x,y,z,rx,1)
    l,m,n = rotatevector(l,m,n,rx,1)
    ux,uy,uz = rotatevector(ux,uy,uz,rx,1)
    
    x -= tx
    y -= ty
    z -= tz
    
    #Update transformation matrices
    if coords is not None:
        #Define rotation and translation matrices
        rotm = rotationM(-rx,-ry,-rz,inverse=True)
        tranm = translationM(tx,ty,tz)
        rotmi = rotationM(-rx,-ry,-rz)
        tranmi = translationM(-tx,-ty,-tz)
        #Dot rotation into forward transform
        coords[0] = np.dot(rotm,coords[0])
        coords[1] = np.dot(np.dot(tranm,rotm),coords[1])
        coords[2] = np.dot(coords[2],rotmi)
        coords[3] = np.dot(coords[3],np.dot(rotmi,tranmi))
    
    return x, y, z, l, m, n, ux, uy, uz

def radgrat(x,y,l,m,n,dpermm,order,wave):
    # Radially grooved grating diffraction
    # Assumes grating in x y plane, with grooves converging at 
    # hubdist in positive y direction
    
    # Don't execute if order or wave is None, there is no way we can reflect
    if order is None or wave is None:
        return x, y, l, m, n
    
    sn = n / np.abs(n)
    d = dpermm * np.sqrt(y**2 + x**2)
    yaw = -np.pi/2 + np.arctan(x/y)
    
    l += np.sin(yaw)*order*wave/d
    m -= np.cos(yaw)*order*wave/d
    with np.errstate(invalid = "ignore"):
        n = sn*np.sqrt(1.0 - l**2 - m**2)
    
    # Check for Evanescence
    with np.errstate(invalid='ignore'):
        removelist = ((l**2 + m**2) > 1)
    
    x[removelist] = np.nan
    y[removelist] = np.nan
    l[removelist] = np.nan
    m[removelist] = np.nan
    n[removelist] = np.nan
    
    return x, y, l, m, n
    
def grat(l,m,n,d,order,wave, eliminate="nan"):
    # Linear grating with groove period d
    # Wavelength wave
    # Groove direction assumed in y direction
    
    # Don't execute if order or wave is None, there is no way we can reflect
    if order is None or wave is None:
        return l, m, n
    
    sn = n / np.abs(n)
    l -= order*wave/d
    with np.errstate(invalid='ignore'):
        n = sn*np.sqrt(1.0 - l**2 - m**2)
    
    # Check for Evanescence
    with np.errstate(invalid='ignore'):
        removelist = ((l**2 + m**2) > 1)

    l[removelist] = np.nan
    m[removelist] = np.nan
    n[removelist] = np.nan
        
    return l, m, n

#Transformation matrix helper functions
def rotationM(rx,ry,rz,inverse=False):
    """Return a rotation matrix, applying rotations in
    X,Y,Z order
    Negate the angle values to be consistent with transform function
    Translation translates the reference frame
    """
    if inverse is True:
        rx,ry,rz = -rx,-ry,-rz
    r1 = tr.rotation_matrix(-rx,[1,0,0])
    r2 = tr.rotation_matrix(-ry,[0,1,0])
    r3 = tr.rotation_matrix(-rz,[0,0,1])
    if inverse is True:
        return np.dot(r1,np.dot(r2,r3))
    else:
        return np.dot(r3,np.dot(r2,r1))

def translationM(tx,ty,tz):
    """
    Return a translation matrix. Negate the values in order
    to be consistent with the transform method.
    Translation translates the reference frame"""
    return tr.translation_matrix([-tx,-ty,-tz])

def applyT(rays,coords,inverse=False):
    """Apply transformation matrix to raylist.
    Only rotations to direction cosines.
    Inverse means going back to global coordinate system.
    Forward means going from global coordinate system to
    local coordinate system.
    """
    i = 0
    if inverse is True:
        i = 2
    #Extract position, wavevector, and surface normals
    on = np.ones(np.shape(rays.x)[0])
    pos = [rays.x,rays.y,rays.z,on]
    wave = [rays.l,rays.m,rays.n,on]
    norm = [rays.ux,rays.uy,rays.uz,on]
    #Apply relevant transformations
    pos = np.dot(coords[i+1],pos)[:3]
    wave = np.dot(coords[i],wave)[:3]
    norm = np.dot(coords[i],norm)[:3]
    x,y,z = pos
    l,m,n = wave
    ux,uy,uz = norm
    #Construct and return new raylist
    return x, y, z, l, m, n, ux, uy, uz 

    
    
    
    
    
