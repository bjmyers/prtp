import numpy as np
import transformationsf as tran
from analyses import analyticYPlane,analyticXPlane,analyticImagePlane

def flat(rays):
    """
    Trace rays to the XY plane
    """
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz
    delta = np.zeros(len(x))
    delta = (-1)*z/n
    ux = np.zeros(len(x))
    uy = np.zeros(len(x))
    uz = np.ones(len(x))
    z = np.zeros(len(x))
    x += delta*l
    y += delta*m
    return x, y, z, l, m, n, ux, uy, uz

def sphere(rays,rad):
    """
    Wrapper for spherical surface.
    """
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz
    
    num = len(x)
    dotol=np.zeros(num)
    mago=np.zeros(num)
    determinant=np.zeros(num)
    
    dotol = l*x + m*y + n*z
    mago = x**2 + y**2 + z**2
    determinant = dotol**2 - mago + rad**2
    
    d1=np.zeros(len(x))
    d2=np.zeros(len(x))
    with np.errstate(invalid='ignore'):
        d1 = np.where(determinant<0,np.nan,(-1)*dotol + np.sqrt(determinant))
        d2 = np.where(determinant<0,np.nan,(-1)*dotol - np.sqrt(determinant))
        d1 = np.where(np.abs(d2) < np.abs(d1),d2,d1)
    
    x = np.where(determinant < 0,np.nan,x+d1*l)
    y = np.where(determinant < 0,np.nan,y+d1*m)
    z = np.where(determinant < 0,np.nan,z+d1*n)
    l = np.where(determinant < 0,np.nan,l)
    m = np.where(determinant < 0,np.nan,m)
    n = np.where(determinant < 0,np.nan,n)
        
    mago = np.sqrt(x**2 + y**2 + z**2)
    with np.errstate(invalid='ignore'):
        ux = np.where(mago==0,np.nan,x/mago)
        uy = np.where(mago==0,np.nan,y/mago)
        uz = np.where(mago==0,np.nan,z/mago)
        
    return x, y, z, l, m, n, ux, uy, uz

def cyl(rays,rad):
    """
    Wrapper for cylindrical surface
    """
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz
    
    num = len(x)
    a=np.zeros(num)
    b=np.zeros(num)
    c=np.zeros(num)
    determinant = np.zeros(num)
    d1=np.zeros(num)
    d2=np.zeros(num)
    mag=np.zeros(num)
    
    a = l**2 + n**2
    b = 2*(x*l+z*n)
    c = x**2 + z**2 - rad**2
    determinant = b**2 - 4*a*c
    
    with np.errstate(invalid='ignore'):
        d1 = np.where(determinant<0,np.nan,(-b + np.sqrt(determinant))/2/a)
        d2 = np.where(determinant<0,np.nan,(-b - np.sqrt(determinant))/2/a)
        d1 = np.where(np.abs(d2) < np.abs(d1),d2,d1)
        
    x = np.where(determinant < 0,np.nan,x+l*d1)
    y = np.where(determinant < 0,np.nan,y+m*d1)
    z = np.where(determinant < 0,np.nan,z+n*d1)
    l = np.where(determinant < 0,np.nan,l)
    m = np.where(determinant < 0,np.nan,m)
    n = np.where(determinant < 0,np.nan,n)
        
    mag = np.sqrt(x**2+z**2)
    with np.errstate(invalid='ignore'):
        ux = np.where(mag==0,0,x/mag)
        uz = np.where(mag==0,0,z/mag)
    uy = np.zeros(num)
        
    return x, y, z, l, m, n, ux, uy, uz


def cylconic(rays,rad,k,maxiter=12):
    """
    Wrapper for cylindrical surface
    """
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz

    delt  = np.ones(len(x))
    delt *= 100
    low   = np.zeros(len(x))
    high  = np.zeros(len(x))
    dL    = np.zeros(len(x))
    dH    = np.zeros(len(x))
    F     = np.zeros(len(x))
    Fx    = np.zeros(len(x))
    Fy    = np.ones(len(x))
    Fp    = np.zeros(len(x))
    
    i = 0
    while (i < maxiter):
        low  = 1 + np.sqrt(np.abs(1 - (1+k)*rad**2*x**2))
        high = rad*x**2
        dL   = (-1)*(1+k)*rad**2*x/np.sqrt(np.abs(1 - (1+k)*rad**2*x**2))
        dH   = 2*rad*x
        F    = y - high/low
        Fx   = (high*dL - low*dH) / (low**2)
        Fp   = Fx*l + Fy*m
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt = (-1)*F/Fp
            x += l*delt
            y += m*delt
            z += n*delt
            if (not (np.abs(delt) > 1.0e-10).any()):
                break
        i += 1
    
    with(np.errstate(invalid='ignore',divide='ignore')):
        unfinishedrays = np.logical_or((np.abs(delt) > 1e-10),(np.isnan(delt)))
    x = np.where(unfinishedrays,np.nan,x)
    y = np.where(unfinishedrays,np.nan,y)
    z = np.where(unfinishedrays,np.nan,z)
    l = np.where(unfinishedrays,np.nan,l)
    m = np.where(unfinishedrays,np.nan,m)
    n = np.where(unfinishedrays,np.nan,n)        
    Fp = np.where(unfinishedrays,np.nan,np.sqrt(Fx**2 + Fy**2))      
    ux = np.where(unfinishedrays,np.nan,Fx/Fp)
    uy = np.where(unfinishedrays,np.nan,Fy/Fp)
    uz = np.where(unfinishedrays,np.nan,0)
        
    return x, y, z, l, m, n, ux, uy, uz


def conic(rays,r,k):
    """Wrapper for conic surface with radius of curvature R
    and conic constant K
    """
    # Temp list so inputs are not modified
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz
    
    num = len(x)
    s = np.zeros(num)
    b = np.zeros(num)
    c = np.zeros(num)
    denom = np.zeros(num)
    disc = np.zeros(num)
    s1 = np.zeros(num)
    s2 = np.zeros(num)
    
    if (k==-1):
        type1rays = (np.abs(n)==1.0)
        s = np.where(n==1,x**2 + y**2 - (2*R*z)/(2*R*n),0)
    else:
        type1rays = [False]*num
        denom = l**2 + m**2 + (k+1)*n**2
        b = (x*l + y*m + ((k+1)*z - r)*n) / denom
        c = (x**2 + y**2 - 2*r*z + (k+1)*z**2) / denom
        disc = b**2 - c
    
    type2rays = (disc > 0)
    noteliminated = np.logical_or(type1rays,type2rays)
    with np.errstate(invalid='ignore'):
        s1 = np.where(noteliminated,-b + np.sqrt(disc),0)
        s2 = np.where(noteliminated,-b - np.sqrt(disc),0)
        s = np.where(noteliminated,np.minimum(np.abs(s1),np.abs(s2)),0)
    
    l = np.where(noteliminated,l,np.nan)
    m = np.where(noteliminated,m,np.nan)
    n = np.where(noteliminated,n,np.nan)
    x = np.where(noteliminated,x + l*s,np.nan)
    y = np.where(noteliminated,y + m*s,np.nan)
    z = np.where(noteliminated,z + n*s,np.nan)
    denom = np.where(noteliminated,np.sqrt(r**2 - k*(x**2+y**2)),1)
    with np.errstate(invalid='ignore'):
        ux = np.where(noteliminated,-x/denom,np.nan)
        uy = np.where(noteliminated,-y/denom,np.nan)
        uz = np.where(noteliminated,-(-r/np.abs(r) * np.sqrt(r**2 - (k+1)*(x**2+y**2)))/denom,np.nan)
        
    return x, y, z, l, m, n, ux, uy, uz
    

def paraxial(rays,f):
    """
    Trace rays through an ideal, paraxial lens.
    Assume optical axis is at xy=0 in z direction
    Surface is in xy plane
    """    
    # Temp list so inputs are not modified
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz
    
    l -= x/f
    m -= y/f
    return x, y, z, l, m, n, ux, uy, uz
    

def paraxialY(rays,f):
    """
    Trace rays through an ideal, paraxial lens.
    Assume optical axis is at xy=0 in z direction
    Surface is in xy plane
    """    
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz
    
    m -= y/f
    return x, y, z, l, m, n, ux, uy, uz
    
    
def torus(rays,rin,rout,eliminate="nan", maxiter=12):
    """Wrapper for toroidal surface. Outer radius
    is in xy plane, inner radius is orthogonal.
    """
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz
    
    num = len(x)
    delt = np.ones(num) * 100
    F = np.zeros(num)
    Fx = np.zeros(num)
    Fy = np.zeros(num)
    Fz = np.zeros(num)
    Fp = np.zeros(num)
    
    i = 0
    while (i < maxiter):
        with(np.errstate(invalid='ignore',divide='ignore')):
            F  = ((z+rin+rout)**2+y**2+x**2+rout**2-rin**2)**2 - (4*rout**2*(y**2+(z+rin+rout)**2))
            Fx =4*x * (-rin**2+(rin+rout+z)**2+rout**2+x**2+y**2)
            Fy = 4*y * (2*rin*(rout+z) + 2*rout*z + z**2+y**2+x**2)
            Fz = 4*(rout+rin+z)*(2*rin*(rout+z) + 2*rout*z+z**2+y**2+x**2)
            Fp = Fx*l + Fy*m + Fz*n
            delt = -F/Fp
            x += l*delt
            y += m*delt
            z += n*delt
            if (not (np.abs(delt) > 1.0e-10).any()):
                break
        i += 1
           
    with(np.errstate(invalid='ignore',divide='ignore')):
        unfinishedrays = (np.abs(delt) > 1e-10) | (np.isnan(delt))
    x = np.where(unfinishedrays,np.nan,x)
    y = np.where(unfinishedrays,np.nan,y)
    z = np.where(unfinishedrays,np.nan,z)
    l = np.where(unfinishedrays,np.nan,l)
    m = np.where(unfinishedrays,np.nan,m)
    n = np.where(unfinishedrays,np.nan,n)
    
    Fp = np.where(unfinishedrays,np.nan,np.sqrt(Fx**2 + Fy**2))
    
    ux = np.where(unfinishedrays,np.nan,Fx/Fp)
    uy = np.where(unfinishedrays,np.nan,Fy/Fp)
    uz = np.where(unfinishedrays,np.nan,Fz/Fp)
       
    return x, y, z, l, m, n, ux, uy, uz


def conicplus(rays,r,k,p,Np,eliminate='nan',maxiter=12):
    """Wrapper for conic surface with radius of curvature R
    and conic constant K
    """
    x,y,z,l,m,n,ux,uy,uz = rays.x,rays.y,rays.z,rays.l,rays.m,rays.n,rays.ux,rays.uy,rays.uz

    num = len(x)
    delt = np.ones(num) * 100
    F = np.zeros(num)
    Fx = np.zeros(num)
    Fy = np.zeros(num)
    Fr = np.zeros(num)
    Fp = np.zeros(num)
    rad = np.zeros(num)
    denom = np.zeros(num)
    a0 = np.zeros(num)
    a1 = np.zeros(num)
    Fz = 1
    c = 1/r
    
    for j in range(1,Np):
        a0 += p[j]*r**(2*j)
        a1 += p[j]*2*j*r**(2*j-1)
    
    i = 0
    while (i < maxiter):
        with(np.errstate(invalid='ignore',divide='ignore')):
            rad = np.sqrt(x**2 + y**2)
            denom = np.sqrt(1-(k+1)*c**2*r**2)+1
            Fr = -(2*c*rad/denom + ((k+1)*rad**3*c**3)/(denom**2*np.sqrt(1-(k+1)*c**2*rad**2)))+a1
            
        F  = z - c*r**2 / denom + a0
        Fx = Fr + x/rad
        Fy = Fr + y/rad
        Fp = Fx*l + Fy*m + Fz*n
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt = -F/Fp
            x += l*delt
            y += m*delt
            z += n*delt
            if (not (np.abs(delt) > 1.0e-10).any()):
                break
        i += 1  

    with(np.errstate(invalid='ignore',divide='ignore')):
        unfinishedrays = (np.abs(delt) > 1e-10) | (np.isnan(delt))
    x = np.where(unfinishedrays,np.nan,x)
    y = np.where(unfinishedrays,np.nan,y)
    z = np.where(unfinishedrays,np.nan,z)
    l = np.where(unfinishedrays,np.nan,l)
    m = np.where(unfinishedrays,np.nan,m)
    n = np.where(unfinishedrays,np.nan,n)
    opd = np.where(unfinishedrays,np.nan,opd)
    Fp = np.where(unfinishedrays,np.nan,np.sqrt(Fx**2 + Fy**2 + Fz**2))
    ux = np.where(unfinishedrays,np.nan,Fx/Fp)
    uy = np.where(unfinishedrays,np.nan,Fy/Fp)
    uz = np.where(unfinishedrays,np.nan,Fz/Fp)

    return x, y, z, l, m, n, ux, uy, uz


def focus(rays,fn,weights=None,nr=None):
    dz1 = fn(rays,weights=weights)
    rays = tran.transform(rays,0,0,dz1,0,0,0)
    rays = flat(rays,nr=nr)
    dz2 = fn(rays,weights=weights)
    rays = tran.transform(rays,0,0,dz2,0,0,0)
    rays = flat(rays,nr=nr)
    
    return (rays,dz1+dz2)

def focusY(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticYPlane,weights=weights,nr=nr)

def focusX(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticXPlane,weights=weights,nr=nr)

def focusI(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticImagePlane,weights=weights,nr=nr)
        
    
