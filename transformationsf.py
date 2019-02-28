import numpy as np
from numba import jit
#'UniTuple(float64[:], 2)(float64[:],float64[:],float64[:],int32,int32)'
@jit(cache=True,nopython=True)
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
    return tempx,tempy,z

#'float[:,:](float64[:,:])'
@jit(cache=True,nopython=True)  
def reflect(rays):
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy

    dot = ux*l + uy*m + uz*n
    l -= 2*dot*ux
    m -= 2*dot*uy
    n -= 2*dot*uz
    
    return [opd, x, y, z, l, m, n, ux, uy, uz]
    
    'UniTuple(float64[:], 2)(float64[:],float64[:],float64[:],int32,int32)'
def refract(rays,n1,n2):
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
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
    
    l,m,n = rotateaxis(l, m, n, cx, cy, cz, t2-t1)
    alpha = np.sqrt(l**2 + m**2 + n**2)
    l = np.where(unfinishedrays,l/alpha,l)
    m = np.where(unfinishedrays,m/alpha,m)
    n = np.where(unfinishedrays,n/alpha,n)
    
    return [opd, x, y, z, l, m, n, ux, uy, uz]
    
@jit(cache=True,nopython=True)
def transform(rays,tx,ty,tz,rx,ry,rz):
    # Coordinate system transform, translations are done first, then rotations in XYZ order
    # Rotations act to rotate a surface via the right hand rule
    # tx,ty,tz are linear translations
    # rx,ry,rz are angular translations
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy

    tx *= -1
    ty *= -1
    tz *= -1
    rx *= -1
    ry *= -1
    rz *= -1    

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
    
    return [opd, x, y, z, l, m, n, ux, uy, uz]
    
    
def itransform(rays,tx,ty,tz,rx,ry,rz):
    # This function inverts the transform function with the same arguments
    # transform(rays,0,10,0,2,0,0) will be undone by itransform(rays,0,10,0,2,0,0)
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    x,y,z = rotatevector(x,y,z,-rz,3)
    l,m,n = rotatevector(l,m,n,-rz,3)
    ux,uy,uz = rotatevector(ux,uy,uz,-rz,3)
    
    x,y,z = rotatevector(x,y,z,-ry,2)
    l,m,n = rotatevector(l,m,n,-ry,2)
    ux,uy,uz = rotatevector(ux,uy,uz,-ry,2)
    
    x,y,z = rotatevector(x,y,z,-rx,1)
    l,m,n = rotatevector(l,m,n,-rx,1)
    ux,uy,uz = rotatevector(ux,uy,uz,-rx,1)
    
    x -= tx
    y -= ty
    z -= tz
    
    return [opd, x, y, z, l, m, n, ux, uy, uz]


def radgrat(rays,wave,dpermm,order,eliminate="nan"):
    # Radially grooved grating diffraction
    # Assumes grating in x y plane, with grooves converging at 
    # hubdist in positive y direction
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    sn = n / np.abs(n)
    d = dpermm * np.sqrt(y**2 + x**2)
    yaw = np.pi/2 - np.arctan(x/np.abs(y))
    
    l += np.sin(yaw)*order*wave/d
    m -= np.cos(yaw)*order*wave/d
    with np.errstate(invalid = "ignore"):
        n = np.sqrt(1.0 - l**2 - m**2)
    
    # Check for Evanescence
    removelist = ((l**2 + m**2) > 1)
    
    if (eliminate.lower() == 'nan'):
        x[removelist] = np.nan
        y[removelist] = np.nan
        z[removelist] = np.nan
        l[removelist] = np.nan
        m[removelist] = np.nan
        n[removelist] = np.nan
        ux[removelist] = np.nan
        uy[removelist] = np.nan
        uz[removelist] = np.nan
        opd[removelist] = np.nan
    
    elif (eliminate.lower() == 'remove'):
        x = x[np.logical_not(removelist)]
        y = y[np.logical_not(removelist)]
        z = z[np.logical_not(removelist)]
        l = l[np.logical_not(removelist)]
        m = m[np.logical_not(removelist)]
        n = n[np.logical_not(removelist)]
        ux = ux[np.logical_not(removelist)]
        uy = uy[np.logical_not(removelist)]
        uz = uz[np.logical_not(removelist)]
        opd = opd[np.logical_not(removelist)]
    
    return [opd, x, y, z, l, m, n, ux, uy, uz]


def radgratW(rays,wave,dpermm,order, eliminate="nan"):
    # Radially grooved grating diffraction with wavelength vector
    # Assumes grating in x y plane, with grooves converging at 
    # hubdist in positive y direction
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    d = dpermm * np.sqrt(y**2 + x**2)
    sn = y/np.abs(y)
    yaw = np.pi/2 + np.arctan(x/np.abs(y))
    
    l += np.sin(yaw)*order*wave/d
    m -= np.cos(yaw)*order*wave/d
    n = np.sqrt(1.0 - l**2 - m**2)
    
    # Check for Evanescence
    removelist = ((l**2 + m**2) > 1)
    
    if (eliminate.lower() == 'nan'):
        x[removelist] = np.nan
        y[removelist] = np.nan
        z[removelist] = np.nan
        l[removelist] = np.nan
        m[removelist] = np.nan
        n[removelist] = np.nan
        ux[removelist] = np.nan
        uy[removelist] = np.nan
        uz[removelist] = np.nan
        opd[removelist] = np.nan
    
    elif (eliminate.lower() == 'remove'):
        x = x[np.logical_not(removelist)]
        y = y[np.logical_not(removelist)]
        z = z[np.logical_not(removelist)]
        l = l[np.logical_not(removelist)]
        m = m[np.logical_not(removelist)]
        n = n[np.logical_not(removelist)]
        ux = ux[np.logical_not(removelist)]
        uy = uy[np.logical_not(removelist)]
        uz = uz[np.logical_not(removelist)]
        opd = opd[np.logical_not(removelist)]
    
    return [opd, x, y, z, l, m, n, ux, uy, uz]
    
    
def grat(rays,d,order,wave, eliminate="nan"):
    # Linear grating with groove period d
    # Wavelength wave
    # Groove direction assumed in y direction
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    l -= order*wave/d
    with np.errstate(invalid='ignore'):
        n = np.sqrt(1.0 - l**2 - m**2)
    
    # Check for Evanescence
    removelist = ((l**2 + m**2) > 1)
    
    if (eliminate.lower() == 'nan'):
        x[removelist] = np.nan
        y[removelist] = np.nan
        z[removelist] = np.nan
        l[removelist] = np.nan
        m[removelist] = np.nan
        n[removelist] = np.nan
        ux[removelist] = np.nan
        uy[removelist] = np.nan
        uz[removelist] = np.nan
        opd[removelist] = np.nan
    
    elif (eliminate.lower() == 'remove'):
        x = x[np.logical_not(removelist)]
        y = y[np.logical_not(removelist)]
        z = z[np.logical_not(removelist)]
        l = l[np.logical_not(removelist)]
        m = m[np.logical_not(removelist)]
        n = n[np.logical_not(removelist)]
        ux = ux[np.logical_not(removelist)]
        uy = uy[np.logical_not(removelist)]
        uz = uz[np.logical_not(removelist)]
        opd = opd[np.logical_not(removelist)]
        
    return [opd, x, y, z, l, m, n, ux, uy, uz]
    
   
    
    
    
    
    
    
