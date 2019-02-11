import numpy as np
import prtp.transformationsf as tran
from prtp.analyses import analyticYPlane,analyticXPlane,analyticImagePlane

def flat(rays,nr=None, eliminate="nan"):
    """
    Trace rays to the XY plane
    """
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    if nr is not None:
        delta = np.zeros(len(x))
        delta = (-1)*z/n
        ux = 0
        uy = 0
        uz = 1
        z = 0
        x += delta*l
        y += delta*m
        opd += delta*nr
    else:
        delta = np.zeros(len(x))
        delta = (-1)*z/n
        ux = 0
        uy = 0
        uz = 1
        z = 0
        x += delta*l
        y += delta*m
    return [opd, x, y, z, l, m, n, ux, uy, uz]

def sphere(rays,rad,nr=None, eliminate="nan"):
    """
    Wrapper for spherical surface.
    """
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    num = len(x)
    dotol=np.zeros(num)
    mago=np.zeros(num)
    determinant=np.zeros(num)
    
    dotol = l*x + m*y + n*z
    mago = x**2 + y**2 + z**2
    determinant = dotol**2 - mago + rad**2
    
    if (eliminate.lower() == "nan"):
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
            
    elif (eliminate.lower() == "remove"):
        removelist = (determinant >= 0)
        x = x[removelist]
        y = y[removelist]
        z = z[removelist]
        l = l[removelist]
        m = m[removelist]
        n = n[removelist]
        ux = ux[removelist]
        uy = uy[removelist]
        uz = uz[removelist]
        dotol = dotol[removelist]
        mago = mago[removelist]
        determinant = determinant[removelist]
        d1 = np.zeros(len(x))
        d2 = np.zeros(len(x))
        d1 = (-1)*dotol + np.sqrt(determinant)
        d2 = (-1)*dotol - np.sqrt(determinant)
        d1 = np.where(np.abs(d2) < np.abs(d1),d2,d1)
        x += d1*l
        y += d1*m
        z += d1*n
        mago = np.sqrt(x**2 + y**2 + z**2)
        with np.errstate(invalid='ignore'):
            ux = np.where(mago==0,np.nan,x/mago)
            uy = np.where(mago==0,np.nan,y/mago)
            uz = np.where(mago==0,np.nan,z/mago)
    
    if nr is not None:
        opd = np.where(determinant < 0,opd,opd+d1*nr)
        
    return [opd, x, y, z, l, m, n, ux, uy, uz]

def cyl(rays,rad,nr=None, eliminate="nan"):
    """
    Wrapper for cylindrical surface
    """
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
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
    
    if (eliminate.lower() == "nan"):
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
        uy = 0
        
    elif (eliminate.lower() == "remove"):
        removelist = (determinant >= 0)
        x = x[removelist]
        y = y[removelist]
        z = z[removelist]
        l = l[removelist]
        m = m[removelist]
        n = n[removelist]
        ux = ux[removelist]
        uy = uy[removelist]
        uz = uz[removelist]
        a = a[removelist]
        b = b[removelist]
        c = c[removelist]
        determinant = determinant[removelist]
        d1=np.zeros(len(x))
        d2=np.zeros(len(x))
        mag=np.zeros(len(x))
        
        d1 = (-b + np.sqrt(determinant))/2/a
        d2 = (-b - np.sqrt(determinant))/2/a
        d1 = np.where(np.abs(d2) < np.abs(d1),d2,d1)
        
        x = x+l*d1
        y = y+m*d1
        z = z+n*d1
        mag = np.sqrt(x**2+z**2)
        with np.errstate(invalid='ignore'):
            ux = np.where(mag==0,0,x/mag)
            uz = np.where(mag==0,0,z/mag)
        uy = 0
    
    if nr is not None:
        opd += d1*nr
        
    return [opd, x, y, z, l, m, n, ux, uy, uz]


def cylconic(rays,rad,k, eliminate="nan", maxiter=12):
    """
    Wrapper for cylindrical surface
    """
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy

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
    
    if (eliminate.lower() == "nan"):
        with(np.errstate(invalid='ignore',divide='ignore')):
            unfinishedrays = (np.abs(delt) > 1e-10) | (np.isnan(delt))
        x = np.where(unfinishedrays,np.nan,x)
        y = np.where(unfinishedrays,np.nan,y)
        z = np.where(unfinishedrays,np.nan,z)
        l = np.where(unfinishedrays,np.nan,l)
        m = np.where(unfinishedrays,np.nan,m)
        n = np.where(unfinishedrays,np.nan,n)
        opd = np.where(unfinishedrays,np.nan,opd)
        
        Fp = np.where(unfinishedrays,np.nan,np.sqrt(Fx**2 + Fy**2))
        
        ux = np.where(unfinishedrays,np.nan,Fx/Fp)
        uy = np.where(unfinishedrays,np.nan,Fy/Fp)
        uz = np.where(unfinishedrays,np.nan,0)
    
    elif (eliminate.lower() == "remove"):
        with(np.errstate(invalid='ignore',divide='ignore')):
            removelist = (np.abs(delt) <= 1e-10) | (np.isnan(delt))
        x = x[removelist]
        y = y[removelist]
        z = z[removelist]
        l = l[removelist]
        m = m[removelist]
        n = n[removelist]
        ux = ux[removelist]
        uy = uy[removelist]
        uz = np.zeros(len(x))
        Fx = Fx[removelist]
        Fy = Fy[removelist]
        Fp = np.zeros(len(x))
        opd = opd[removelist]
        
        Fp = np.sqrt(Fx**2 * Fy**2)
        ux = Fx / Fp
        uy = Fy / Fp
        
    return [opd, x, y, z, l, m, n, ux, uy, uz]


def conic(rays,r,k,nr=None, eliminate="nan"):
    """Wrapper for conic surface with radius of curvature R
    and conic constant K
    """
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
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
    noteliminated = (type1rays | type2rays)
    with np.errstate(invalid='ignore'):
        s1 = np.where(noteliminated,-b + np.sqrt(disc),0)
        s2 = np.where(noteliminated,-b - np.sqrt(disc),0)
        s = np.where(noteliminated,np.minimum(np.abs(s1),np.abs(s2)),0)
    
    if (eliminate.lower() == "nan"):
        l = np.where(noteliminated,l,np.nan)
        m = np.where(noteliminated,m,np.nan)
        n = np.where(noteliminated,n,np.nan)
        x = np.where(noteliminated,x + l*s,np.nan)
        y = np.where(noteliminated,y + m*s,np.nan)
        z = np.where(noteliminated,z + n*s,np.nan)
        if (nr is not None):
            opd = np.where(noteliminated,opd + s*nr,np.nan)
        denom = np.where(noteliminated,np.sqrt(r**2 - k*(x**2+y**2)),1)
        with np.errstate(invalid='ignore'):
            ux = np.where(noteliminated,-x/denom,np.nan)
            uy = np.where(noteliminated,-y/denom,np.nan)
            uz = np.where(noteliminated,-(-r/np.abs(r) * np.sqrt(r**2 - (k+1)*(x**2+y**2)))/denom,np.nan)
        
    elif (eliminate.lower() == "remove"):
        removelist = (s!=0)
        x = x[removelist]
        y = y[removelist]
        z = z[removelist]
        l = l[removelist]
        m = m[removelist]
        n = n[removelist]
        ux = ux[removelist]
        uy = uy[removelist]
        uz = uz[removelist]
        denom = denom[removelist]
        s = s[removelist]
        opd = opd[removelist]
        x += l*s
        y += m*s
        z += n*s
        if (nr is not None):
            opd += s*nr
        ux = -x/denom
        uy = -y/denom
        with np.errstate(invalid='ignore'):
            uz = -(-r/np.abs(r) * np.sqrt(r**2 - (k+1)*(x**2+y**2)))/denom
        
    return [opd, x, y, z, l, m, n, ux, uy, uz]
    

def paraxial(rays,f,eliminate="nan"):
    """
    Trace rays through an ideal, paraxial lens.
    Assume optical axis is at xy=0 in z direction
    Surface is in xy plane
    """    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    l -= x/f
    m -= y/f
    return [opd, x, y, z, l, m, n, ux, uy, uz]
    

def paraxialY(rays,f,eliminate="nan"):
    """
    Trace rays through an ideal, paraxial lens.
    Assume optical axis is at xy=0 in z direction
    Surface is in xy plane
    """    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    m -= y/f
    return [opd, x, y, z, l, m, n, ux, uy, uz]
    
    
def torus(rays,rin,rout,eliminate="nan", maxiter=12):
    """Wrapper for toroidal surface. Outer radius
    is in xy plane, inner radius is orthogonal.
    """
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
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
           
    if (eliminate.lower() == "nan"):
        with(np.errstate(invalid='ignore',divide='ignore')):
            unfinishedrays = (np.abs(delt) > 1e-10) | (np.isnan(delt))
        x = np.where(unfinishedrays,np.nan,x)
        y = np.where(unfinishedrays,np.nan,y)
        z = np.where(unfinishedrays,np.nan,z)
        l = np.where(unfinishedrays,np.nan,l)
        m = np.where(unfinishedrays,np.nan,m)
        n = np.where(unfinishedrays,np.nan,n)
        opd = np.where(unfinishedrays,np.nan,opd)
        
        Fp = np.where(unfinishedrays,np.nan,np.sqrt(Fx**2 + Fy**2))
        
        ux = np.where(unfinishedrays,np.nan,Fx/Fp)
        uy = np.where(unfinishedrays,np.nan,Fy/Fp)
        uz = np.where(unfinishedrays,np.nan,Fz/Fp)
        
    elif (eliminate.lower() == "remove"):
        with(np.errstate(invalid='ignore',divide='ignore')):
            removelist = (np.abs(delt) <= 1e-10)
        x = x[removelist]
        y = y[removelist]
        z = z[removelist]
        l = l[removelist]
        m = m[removelist]
        n = n[removelist]
        ux = ux[removelist]
        uy = uy[removelist]
        uz = uz[removelist]
        Fx = Fx[removelist]
        Fy = Fy[removelist]
        Fz = Fz[removelist]
        Fp = np.zeros(len(x))
        opd = opd[removelist]
        
        Fp = np.sqrt(Fx**2 * Fy**2)
        ux = Fx / Fp
        uy = Fy / Fp
        uz = Fz / Fp
        
    return [opd, x, y, z, l, m, n, ux, uy, uz]


def conicplus(rays,r,k,p,Np,nr=None,eliminate='nan',maxiter=12):
    """Wrapper for conic surface with radius of curvature R
    and conic constant K
    """
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
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

    if (eliminate.lower() == "nan"):
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
        if (nr is not None):
            opd = np.where(unfinishedrays,np.nan,opd + nr*delt)

    elif (eliminate.lower() == "remove"):
        with(np.errstate(invalid='ignore',divide='ignore')):
            removelist = (np.abs(delt) <= 1e-10)
        x = x[removelist]
        y = y[removelist]
        z = z[removelist]
        l = l[removelist]
        m = m[removelist]
        n = n[removelist]
        ux = ux[removelist]
        uy = uy[removelist]
        uz = uz[removelist]
        Fx = Fx[removelist]
        Fy = Fy[removelist]
        Fp = np.zeros(len(x))
        opd = opd[removelist]
        delt = delt[removelist]
        if (nr is not None):
            opd += nr*delt
        
        Fp = np.sqrt(Fx**2 + Fy**2 + Fz**2)
        ux = Fx/Fp
        uy = Fy/Fp
        uz = Fz/Fp

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def focus(rays,fn,weights=None,nr=None):
    dz1 = fn(rays,weights=weights)
    rays = tran.transform(rays,0,0,-dz1,0,0,0)
    rays = flat(rays,nr=nr)
    dz2 = fn(rays,weights=weights)
    rays = tran.transform(rays,0,0,-dz2,0,0,0)
    rays = flat(rays,nr=nr)
    
    return (rays,dz1+dz2)

def focusY(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticYPlane,weights=weights,nr=nr)

def focusX(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticXPlane,weights=weights,nr=nr)

def focusI(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticImagePlane,weights=weights,nr=nr)
