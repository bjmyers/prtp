import numpy as np
import prtp.specialFunctions as funcs

def wolterprimary(rays,r0,z0,psi,nr=None,eliminate="nan",maxiter=10):
    """
    This function traces to a Wolter I primary mirror
    Defined by Van Speybroeck prescription
    For WFS test, use flat to get rays close so they find correct intersection
    Surface should be placed at common focus with z+ pointing toward mirrors
    """
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    # Van Speybroeck Parameters
    alpha = 0.25*np.arctan(r0/z0)
    thetah = 2*(1+2*psi)/(1+psi)*alpha
    thetap = 2*psi/(1+psi) * alpha
    p = z0*np.tan(4*alpha)*np.tan(thetap)
    d = z0*np.tan(4*alpha)*np.tan(4*alpha-thetah)
    e = np.cos(4*alpha)*(1+np.tan(4*alpha)*np.tan(thetah))
    
    # Trace rays to mirror
    delt = np.ones(len(x))*100
    Fz = 2*p
    i = 0
    while (i < maxiter):
        F = 2*p*z + p**2 + 4*e**2*p*d/(e**2-1) - x[:]**2 - y[:]**2
        Fx = -2.*x
        Fy = -2.*y
        Fp = Fx*l + Fy*m + Fz*n
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt = (-1)*F/Fp
            x += l*delt
            y += m*delt
            z += n*delt
            if nr is not None:
                opd += nr*delt
            if (not (np.abs(delt) > 1.0e-10).any()):
                break
        i += 1
    
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx**2 + Fy**2 + Fz**2)
        ux = Fx/Fp
        uy = Fy/Fp
        uz = Fz/Fp

    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-9,eliminate = eliminate)
    
    
def woltersecondary(rays,r0,z0,psi,eliminate="nan",maxiter=10):
    """
    This function traces to a Wolter I secondary mirror
    Defined by Van Speybroeck prescription
    For WFS test, use flat to get rays close so they find correct intersection
    Surface should be placed at common focus with z+ pointing toward mirrors
    """

    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    # Van Speybroeck Parameters
    alpha = 0.25*np.arctan(r0/z0)
    thetah = 2*(1+2*psi)/(1+psi)*alpha
    thetap = 2*psi/(1+psi) * alpha
    p = z0*np.tan(4*alpha)*np.tan(thetap)
    d = z0*np.tan(4*alpha)*np.tan(4*alpha-thetah)
    e = np.cos(4*alpha)*(1+np.tan(4*alpha)*np.tan(thetah))
    
    # Trace rays to mirror
    delt = np.ones(len(x))*100
    i = 0
    while (i < maxiter):
        F = e**2*(d+z)**2 - z**2 - x**2 - y**2
        Fx = -2.*x
        Fy = -2.*y
        Fz = 2*e**2*(d+z) - 2*z
        Fp = Fx*l + Fy*m + Fz*n 
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F/Fp
            x += l*delt
            y += m*delt
            z += n*delt
            if (not (np.abs(delt) > 1.0e-10).any()):
                break
        i += 1
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx**2 + Fy**2 + Fz**2)
        ux = Fx/Fp
        uy = Fy/Fp
        uz = Fz/Fp
    return [opd,x,y,z,l,m,n,ux,uy,uz]
    #return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-9,eliminate = eliminate)
    
    
def woltersine(rays,r0,z0,amp,freq,eliminate="nan",maxiter=12):
    """
    This function traces to a Wolter I primary mirror with sinusoidal perturbation
    Defined by Van Speybroeck prescription
    For WFS test, use flat to get rays close so they find correct intersection
    Surface should be placed at common focus with z+ pointing toward mirrors
    """
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    alpha = .25*np.arctan(r0/z0)
    thetah = 3.*alpha
    thetap = alpha
    p = z0*np.tan(4*alpha)*np.tan(thetap)
    d = z0*np.tan(4*alpha)*np.tan(4*alpha-thetah)
    e = np.cos(4*alpha)*(1+np.tan(4*alpha)*np.tan(thetah))
    
    # Trace rays to mirror
    delt = np.ones(len(x))*100
    i = 0
    while (i < maxiter):
        rad =  np.sqrt(x[:]**2+y[:]**2) + amp*np.sin(2*np.pi*freq*z[:])
        F = 2*p*z[:] + p**2 + 4*e**2*p*d/(e**2-1) - rad**2
        Fx = -2.*x[:]
        Fy = -2.*y[:]
        Fz = 2.*p - 2*rad*amp*2*np.pi*freq*np.cos(2*np.pi*freq*z[:])
        Fp = Fx[:]*l[:] + Fy[:]*m[:] + Fz[:]*n[:]
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F[:]/Fp[:]
            x[:] += l[:]*delt[:]
            y[:] += m[:]*delt[:]
            z[:] += n[:]*delt[:]
            if (not (np.abs(delt) > 1.0e-10).any()):
                break
        i += 1
    
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx[:]**2 + Fy[:]**2 + Fz[:]**2)
        ux[:] = Fx[:]/Fp[:]
        uy[:] = Fy[:]/Fp[:]
        uz[:] = Fz[:]/Fp[:]

    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-9,eliminate = eliminate)
    

def wolterprimLL(rays,r0,z0,zmax,zmin,dphi,coeff,axial,az,cnum,eliminate='nan',maxiter=10):
    """
    Construct a paraboloid as in Wolter but with Legendre-Legendre
    deformations. Define Legendre and Legendre derivative functions.
    Pass in coeff, axial order, and azimuthal order as in Zemax implementation
    """
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    alpha = .25*np.arctan(r0/z0)
    thetah = 3.*alpha
    thetap = alpha
    p = z0*np.tan(4*alpha)*np.tan(thetap)
    d = z0*np.tan(4*alpha)*np.tan(4*alpha-thetah)
    e = np.cos(4*alpha)*(1+np.tan(4*alpha)*np.tan(thetah))
    
    delt = np.ones(len(x))*100
    i = 0
    num = len(x)
    add = np.zeros(num)
    addx = np.zeros(num)
    addy = np.zeros(num)
    addz = np.zeros(num)
    
    while (i < maxiter):
        
        ang = np.arctan2(y,x)
        targ = 2*ang/dphi
        zarg = (z - ((zmax+zmin)/2.)) / ((zmax-zmin)/2.)
        add *= 0
        addx *= 0
        addy *= 0 
        addz *= 0
        
        legre1 = funcs.legendre(zarg,axial)
        legre2 = funcs.legendre(targ,az)
        legrep1 = funcs.legendrep(targ,az)
        legrep2 = funcs.legendrep(zarg,axial)
        for a in range(cnum): #This for loop can probably be condensed using numpy
            add += coeff[a] * legre1 * legre2
            addx -= coeff[a] * legre1 * legrep1 * (2/dphi)*(y[:]/(y[:]**2+x[:]**2))
            addy += coeff[a] * legre1 * legrep1 * (2/dphi)*(y[:]/(y[:]**2+x[:]**2))
            addz += coeff[a] * legrep2 * legre2 * 2 / (zmax - zmin)
        
        G = np.sqrt(x**2 + y**2) + add
        F = -(G**2 - p**2 - 2*p*z - 4*e**2*p*d/(e**2-1))
        Fx = -2*G*(x/np.sqrt(x**2 + y**2) + addy)
        Fy = -2*G*(y/np.sqrt(x**2 + y**2) + addy)
        Fz = 2*p - 2*G*addz
        Fp = Fx*l + Fy*m + Fz*n
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F[:]/Fp[:]
            x[:] += l[:]*delt[:]
            y[:] += m[:]*delt[:]
            z[:] += n[:]*delt[:]
            if (not (np.abs(delt) > 1.0e-10).any()):
                break
        i += 1
    
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx[:]**2 + Fy[:]**2 + Fz[:]**2)
        ux[:] = Fx[:]/Fp[:]
        uy[:] = Fy[:]/Fp[:]
        uz[:] = Fz[:]/Fp[:]
    
    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-9,eliminate = eliminate)


def woltersecLL(rays,r0,z0,psi,zmax,zmin,dphi,coeff,axial,az,cnum,eliminate='nan',maxiter=10):
    """
    Construct a hyperboloid as in Wolter but with Legendre-Legendre
    deformations. Define Legendre and Legendre derivative functions.
    Pass in coeff, axial order, and azimuthal order as in Zemax implementation
    """
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    alpha = .25*np.arctan(r0/z0)
    thetah = 3.*alpha
    thetap = alpha
    p = z0*np.tan(4*alpha)*np.tan(thetap)
    d = z0*np.tan(4*alpha)*np.tan(4*alpha-thetah)
    e = np.cos(4*alpha)*(1+np.tan(4*alpha)*np.tan(thetah))
    
    num = len(x)
    delt = np.ones(num)*100
    i = 0
    add = np.zeros(num)
    addx = np.zeros(num)
    addy = np.zeros(num)
    addz = np.zeros(num)
    
    while (i < maxiter):
        
        ang = np.arctan2(y,x)
        zarg = (z - ((zmax + zmin)/2.0)) / ((zmax-zmin)/2.0)
        targ = 2*ang/dphi
        add *= 0
        addx *= 0
        addy *= 0
        addz *= 0
        
        legre1 = funcs.legendre(zarg,axial)
        legre2 = funcs.legendre(targ,az)
        legrep1 = funcs.legendrep(targ,az)
        legrep2 = funcs.legendrep(zarg,axial)
        for a in range(cnum): #This for loop can probably be condensed using numpy
            add += coeff[a] * legre1 * legre2
            addx -= coeff[a] * legre1 * legrep1 * (2/dphi)*(y[:]/(y[:]**2+x[:]**2))
            addy += coeff[a] * legre1 * legrep1 * (2/dphi)*(y[:]/(y[:]**2+x[:]**2))
            addz += coeff[a] * legrep2 * legre2 * 2 / (zmax - zmin)
        
        G = np.sqrt(x**2 + y**2) + add
        F = -(G**2 - e**2*(d+z)**2 + z**2)
        Fx = -2*G*(x/np.sqrt(x**2+y**2)+addx)
        Fy = -2*G*(y/np.sqrt(x**2+y**2)+addy)
        Fz = 2*e**2*(d+z) - 2*z - 2*G*addz
        Fp = Fx*l + Fy*m + Fz*n
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F[:]/Fp[:]
            x[:] += l[:]*delt[:]
            y[:] += m[:]*delt[:]
            z[:] += n[:]*delt[:]
            if (not (np.abs(delt) > 1.0e-7).any()):
                break
        i += 1
    
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx[:]**2 + Fy[:]**2 + Fz[:]**2)
        ux[:] = Fx[:]/Fp[:]
        uy[:] = Fy[:]/Fp[:]
        uz[:] = Fz[:]/Fp[:]
    
    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-7,eliminate = eliminate)


def wsprimary(rays,alpha,z0,psi,eliminate='nan',maxiter=10):
    '''
    This function traces to a Wolter-Schwarzschild primary mirror
    Defined by Van Speybroeck prescription
    Surface should be placed at common focus with z+ pointing toward mirrors
    If ray is within inner radius of mirror (defined by betas), it will be
    traced to minimum z position
    '''
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    # Chase parameters
    betas = 4*alpha
    ff = z0/np.cos(betas)
    g = ff/psi
    k = np.tan(betas/2)**2
    
    num = len(x)
    delt = np.ones(num)*100
    kterm = np.zeros(num)
    r = np.zeros(num)
    i = 0
    
    while (i < maxiter):
        
        x2y2 = x**2 + y**2
        
        beta = np.arcsin(np.sqrt(x2y2)/ff)
        
        case1 = (beta <= betas)
        beta = np.where(case1,betas,beta)
        kterm = np.where(np.logical_not(case1),(1/k)*np.tan(beta/2)**2-1,0)
        
        F = -z - ff*np.sin(betas/2)**2 + ff**2*np.sin(beta)**2/(4*ff*np.sin(betas/2)**2) + g*np.cos(beta/2)**4*(kterm)**(1-k)
        Fb = ff**2*np.sin(beta)*np.cos(beta)/(2*ff*np.sin(betas/2)**2) - 2*g*np.cos(beta/2)**3*np.sin(beta/2)*(kterm)**(1-k) + g*(1-k)*np.cos(beta/2)*np.sin(beta/2)*(kterm)**(-k)*(1/k)
        Fz[:] = -1
        
        r = np.where(case1,np.sqrt(x2y2),0)
        Fb = np.where(case1,ff*2*np.sin(betas)*np.cos(betas)/(2*ff*np.sin(betas/2)**2) + g*(1-k)*np.cos(betas/2)*np.sin(betas/2)*(1-k),Fb[:])
        F = np.where(case1,F + (r - ff*npsin(betas))*z/(r**2 + z**2)**2*Fb,F)
        Fz = np.where(case1,Fz + (r-ff*np.sin(betas))*(r**2-z**2)/(r**2 + z**2)**2*Fb,Fz)
        
        dbdx = x/np.sqrt(1-(x2y2)/ff**2)/ff/np.sqrt(x2y2)
        dbdy = y/np.sqrt(1-(x2y2)/ff**2)/ff/np.sqrt(x2y2)
        Fx = Fb * dbdx
        Fy = Fb * dbdy
        Fp = Fx*l + Fy*m - n
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F[:]/Fp[:]
            x[:] += l[:]*delt[:]
            y[:] += m[:]*delt[:]
            z[:] += n[:]*delt[:]
            if (not (np.abs(delt) > 1.0e-8).any()):
                break
        i += 1
    
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx[:]**2 + Fy[:]**2 + Fz[:]**2)
        ux[:] = -Fx[:]/Fp[:]
        uy[:] = -Fy[:]/Fp[:]
        uz[:] = -Fz[:]/Fp[:]
    
    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-8,eliminate = eliminate)


def wssecondary(rays,alpha,z0,psi,eliminate='nan',maxiter=10):
    '''
    This function traces to a Wolter-Schwarzschild secondary mirror
    Defined by Van Speybroeck prescription
    Surface should be placed at common focus with z+ pointing toward mirrors
    If ray is within inner radius of mirror (defined by betas), it will be
    traced to minimum z position
    '''
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    # Chase parameters
    betas = 4*alpha
    ff = z0/np.cos(betas)
    g = ff/psi
    k = np.tan(betas/2)**2
    
    num = len(x)
    delt = np.ones(num)*100
    kterm = np.zeros(num)
    a = np.zeros(num)
    r = np.zeros(num)
    i = 0
    
    while (i < maxiter):
        
        x2y2 = x**2 + y**2
        
        beta = np.atan2(np.sqrt(x2y2),z)
        
        case1 = (beta <= betas)
        beta = np.where(case1,betas,beta)
        a = np.where(case1,1/ff,(1-np.cos(beta))/(1-np.cos(betas))/ff + (1+np.cos(beta))/(2*g)*(kterm)**(1+k))
        kterm = np.where(case1,0,(1/k)*np.tan(beta/2)**2-1)
        
        F = -z + np.cos(beta)/a
        
        dadb = np.where(case1,0,np.sin(beta)/ff/(1-np.cos(betas)) - np.sin(beta)/(2*g)*(kterm)**(1+k) + (k+1)*(np.cos(beta) + 1) * np.tan(beta/2)*kterm**k/2/g/k/(np.cos(beta/2)**2))
        Fb = np.where(case1,0,-np.sin(beta)/a - np.cos(beta)/a**2*dadb)
        dadbs = np.where(case1,np.sin(betas)/ff/(1-np.cos(betas)) + (k+1)*(np.cos(betas)+1)*np.tan(betas/2)/np.cos(betas/2)**2/2/g/k,0)
        dbdzs = np.where(case1, -1*np.sin(betas)**2/np.sqrt(x2y2),0)
        gam = np.where(case1,(-ff*np.sin(betas) - ff**2*np.cos(betas)*dadbs)*dbdzs,0)
        dbdx = np.where(case1,0,x*z/(x2y2 + z**2)/np.sqrt(x2y2))
        dbdy = np.where(case1,0,y*z/(x2y2 + z**2)/np.sqrt(x2y2))
        dbdz = np.where(case1,0,-np.sqrt(x2y2)/(x2y2 + z**2))
        
        F = np.where(case1,F + gam*(z - np.sqrt(x2y2)/np.tan(betas)),F)
        Fx = np.where(case1,-2.0/np.tan(betas)*x/np.sqrt(x2y2),Fb * dbdx)
        Fx = np.where(case1,-2.0/np.tan(betas)*y/np.sqrt(x2y2),Fb * dbdy)
        Fz = np.where(case1,gam - 1,-1.0 + Fb*dbdz)
        
        Fp = Fx*l + Fy*m + Fz*n
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F[:]/Fp[:]
            x[:] += l[:]*delt[:]
            y[:] += m[:]*delt[:]
            z[:] += n[:]*delt[:]
            if (not (np.abs(delt) > 1.0e-8).any()):
                break
    
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx[:]**2 + Fy[:]**2 + Fz[:]**2)
        ux[:] = Fx[:]/Fp[:]
        uy[:] = Fy[:]/Fp[:]
        uz[:] = Fz[:]/Fp[:]
    
    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-8,eliminate = eliminate)


def spoCone(rays,R0,tg, eliminate='nan'):
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    sl = np.tan(tg)
    
    # Quadtratic Formula for Ray Advancement Distance
    A = n**2*sl**2 - m**2 - l**2
    B = 2*n*sl*R0 + 2*z*sl**2*n - 2*x*l - 2*y*m
    C = R0**2 + 2*sl*R0*z + z**2*sl**2 - x**2 - y**2
    det = B**2 - 4*A*C
    
    case1 = (det >= 0)
    with np.errstate(invalid='ignore'):
        t1 = np.where(case1,(-B + np.sqrt(det))/(2*A),0)
        t2 = np.where(case1,(-B - np.sqrt(det))/(2*A),0)
        t1 = np.where(np.abs(t2) < np.abs(t1),t2,t1)
        x += np.where(case1,t1*l,x)
        y += np.where(case1,t1*m,y)
        z += np.where(case1,t1*nmz)
        ux = np.where(case1,-x/np.sqrt(x**2 + y**2)*np.cos(tg),ux)
        uy = np.where(case1,-y/np.sqrt(x**2 + y**2)*np.cos(tg),uy)
        uz = np.where(case1,np.sin(tg),uz)
    
    # Determine a boolean array of rays which missed
    missed = np.ones(len(x))
    missed[case1] = 0
    
    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],missed,thresh = 1e-2,eliminate = eliminate)


def ellipsoidWoltLL(rays,r0,z0,psi,S,zmax,zmin,dphi,coeff,axial,ax,cnum,eliminate='nan',maxiter=10):
    '''
    Construct an ellipsoid primary but with Legendre-Legendre
    deformations. Define Legendre and Legendre derivative functions.
    Pass in coeff, axial order, and azimuthal order as in Zemax implementation
    '''
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    # Telescope Parameters
    P = R0/np.sin((psi*np.arcsin(R0/z0) - np.arcsin(R0/S))/(1+psi))
    ff = (S + P)/2
    bq = -(R0**2 + (ff-P)**2 + ff**2)
    cq = ff**2*(ff-P)**2
    aa = np.sqrt((-bq + np.sqrt(bq**2-4*cq))/2.0)
    bb = np.sqrt(aa**2 - ff**2)
    e = ff/aa
    zfoc = ff-P+z0
    
    num = len(x)
    delt = np.ones(num)*100
    i = 0
    add = np.zeros(num)
    addx = np.zeros(num)
    addy = np.zeros(num)
    addz = np.zeros(num)
    
    while (i < maxiter):
        
        ang = np.arctan2(y,x)
        zarg = (z - ((zmax + zmin)/2.0)) / ((zmax-zmin)/2.0)
        targ = 2*ang/dphi
        add *= 0
        addx *= 0
        addy *= 0
        addz *= 0
        
        legre1 = funcs.legendre(zarg,axial)
        legre2 = funcs.legendre(targ,az)
        legrep1 = funcs.legendrep(targ,az)
        legrep2 = funcs.legendrep(zarg,axial)
        for a in range(cnum): #This for loop can probably be condensed using numpy
            add += coeff[a] * legre1 * legre2
            addx -= coeff[a] * legre1 * legrep1 * (2/dphi)*(y[:]/(y[:]**2+x[:]**2))
            addy += coeff[a] * legre1 * legrep1 * (2/dphi)*(y[:]/(y[:]**2+x[:]**2))
            addz += coeff[a] * legrep2 * legre2 * 2 / (zmax - zmin)
        
        G = np.sqrt(x**2 + y**2) + add
        F = (z-zfoc)**2/aa**2 + G**2/bb**2 - 1
        Fx = 2*G/bb**2*(x/np.sqrt(x**2 + y**2) + addx)
        Fy = 2*G/bb**2*(y/np.sqrt(x**2 + y**2) + addy)
        Fz = 2*(z-zfoc)/aa**2 + (2*G/bb**2)*addz
        Fp = Fx*l + Fy*m + Fz*n
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F[:]/Fp[:]
            x[:] += l[:]*delt[:]
            y[:] += m[:]*delt[:]
            z[:] += n[:]*delt[:]
            if (not (np.abs(delt) > 1.0e-10).any()):
                break
    
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx[:]**2 + Fy[:]**2 + Fz[:]**2)
        ux[:] = Fx[:]/Fp[:]
        uy[:] = Fy[:]/Fp[:]
        uz[:] = Fz[:]/Fp[:]
    
    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-9,eliminate = eliminate)


def wsprimaryBack(rays,alpha,z0,psi,thick):
    '''
    This function traces to a Wolter-Schwarzschild primary mirror
    Defined by Van Speybroeck prescription
    Surface should be placed at common focus with z+ pointing toward mirrors
    If ray is within inner radius of mirror (defined by betas), it will be
    traced to minimum z position
    '''
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    # Chase parameters
    betas = 4*alpha
    ff = z0/np.cos(betas)
    g = ff/psi
    k = np.tan(betas/2)**2
    
    num = len(x)
    delt = np.ones(num)*100
    kterm = np.zeros(num)
    a = np.zeros(num)
    r = np.zeros(num)
    i = 0
    
    while (i < maxiter):
        
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y,x)
        x2 = (r-thick)*np.cos(theta)
        y2 = (r-thick)*np.sin(theta)
        
        beta = np.arcsin(np.sqrt(x2**2 + y2**2)/ff)
        case1 = (beta <= betas)
        beta = np.where(case1,betas,beta)
        kterm = np.where(case1,0,(1/k)*np.tan(beta/2)**2 - 1)
        
        F = -z - ff*np.sin(betas/2)**2 + ff**2*np.sin(beta)**2/(4*ff*np.sin(betas/2)**2) + g*np.cos(beta/2)**4*(kterm)**(1-k)
        Fb = ff**2*np.sin(beta)*np.cos(beta)/(2*ff*np.sin(betas/2)**2) - 2*g*np.cos(beta/2)**3*np.sin(beta/2)*(kterm)**(1-k) + g*(1-k)*np.cos(beta/2)*np.sin(beta/2)*(kterm)**(-k)*(1/k)
        
        r = np.where(case1,np.sqrt(x2**2 + y2**2),0)
        Fb = np.where(case1,ff**2*mp.sin(betas)*np.cos(betas)/(2*ff*np.sin(betas/2)**2) + g*(1-k)*np.cos(betas/2)*np.sin(betas/2)*(1/k),Fb)
        F = np.where(case1,F + (r - ff*np.sin(betas))*z(i)/(r**2+z(i)**2)*Fb,F)
        Fz = np.where(case1,Fz + (r-ff*np.sin(betas))*(r**2-z(i)**2)/(r**2+z(i)**2)**2*Fb,-1)
        
        dbdx = x2/np.sqrt(1-(x2**2+y2**2)/ff**2)/ff/np.sqrt(x2**2+y2**2)
        dbdy = y2/np.sqrt(1-(x2**2+y2**2)/ff**2)/ff/np.sqrt(x2**2+y2**2)
        Fx = Fb * dbdx
        Fy = Fb * dbdy
        Fp = Fx*l + Fy*m + Fz*n
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F[:]/Fp[:]
            x[:] += l[:]*delt[:]
            y[:] += m[:]*delt[:]
            z[:] += n[:]*delt[:]
            if (not (np.abs(delt) > 1.0e-8).any()):
                break
                
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx[:]**2 + Fy[:]**2 + Fz[:]**2)
        ux[:] = -Fx[:]/Fp[:]
        uy[:] = -Fy[:]/Fp[:]
        uz[:] = -Fz[:]/Fp[:]
    
    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-8,eliminate = eliminate)


def wssecondaryBack(rays,alpha,z0,psi,thick,eliminate='nan',maxiter=10):
    '''
    This function traces to a Wolter-Schwarzschild secondary mirror
    Defined by Van Speybroeck prescription
    Surface should be placed at common focus with z+ pointing toward mirrors
    If ray is within inner radius of mirror (defined by betas), it will be
    traced to minimum z position
    '''
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    # Chase parameters
    betas = 4*alpha
    ff = z0/np.cos(betas)
    g = ff/psi
    k = np.tan(betas/2)**2
    
    num = len(x)
    delt = np.ones(num)*100
    kterm = np.zeros(num)
    a = np.zeros(num)
    r = np.zeros(num)
    i = 0
    
    while (i < maxiter):
        
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y,x)
        x2 = (r-thick)*np.cos(theta)
        y2 = (r-thick)*np.sin(theta)
        x2y2 = x2**2 + y2**2
        
        beta = np.arctan2(np.sqrt(x2**2 + y2**2),z)
        case1 = (beta <= betas)
        beta = np.where(case1,betas,beta)
        kterm = np.where(case1,0,(1/k)*np.tan(beta/2)**2 - 1)
        a = np.where(case1,1/ff,(1-np.cos(beta))/(1-np.cos(betas))/ff + (1+np.cos(beta))/(2*g)*(kterm)**(1+k))
        
        F = -z + np.cos(beta)/a
        
        dadb = np.where(case1,0,np.sin(beta)/ff/(1-np.cos(betas)) - np.sin(beta)/(2*g)*(kterm)**(1+k) + (k+1)*(np.cos(beta) + 1) * np.tan(beta/2)*kterm**k/2/g/k/(np.cos(beta/2)**2))
        Fb = np.where(case1,0,-np.sin(beta)/a - np.cos(beta)/a**2*dadb)
        dadbs = np.where(case1,np.sin(betas)/ff/(1-np.cos(betas)) + (k+1)*(np.cos(betas)+1)*np.tan(betas/2)/np.cos(betas/2)**2/2/g/k,0)
        dbdzs = np.where(case1, -1*np.sin(betas)**2/np.sqrt(x2y2),0)
        gam = np.where(case1,(-ff*np.sin(betas) - ff**2*np.cos(betas)*dadbs)*dbdzs,0)
        dbdx = np.where(case1,0,x2*z/(x2y2 + z**2)/np.sqrt(x2y2))
        dbdy = np.where(case1,0,y2*z/(x2y2 + z**2)/np.sqrt(x2y2))
        dbdz = np.where(case1,0,-np.sqrt(x2y2)/(x2y2 + z**2))
        
        F = np.where(case1,F + gam*(z - np.sqrt(x2y2)/np.tan(betas)),F)
        Fx = np.where(case1,-2.0/np.tan(betas)*x2/np.sqrt(x2y2),Fb * dbdx)
        Fx = np.where(case1,-2.0/np.tan(betas)*y2/np.sqrt(x2y2),Fb * dbdy)
        Fz = np.where(case1,gam - 1,-1.0 + Fb*dbdz)
        
        Fp = Fx*l + Fy*m + Fz*n
        
        with(np.errstate(invalid='ignore',divide='ignore')):
            delt[:] = (-1)*F[:]/Fp[:]
            x[:] += l[:]*delt[:]
            y[:] += m[:]*delt[:]
            z[:] += n[:]*delt[:]
            if (not (np.abs(delt) > 1.0e-8).any()):
                break
    
    with(np.errstate(invalid='ignore')):
        Fp = np.sqrt(Fx[:]**2 + Fy[:]**2 + Fz[:]**2)
        ux[:] = Fx[:]/Fp[:]
        uy[:] = Fy[:]/Fp[:]
        uz[:] = Fz[:]/Fp[:]
    
    return funcs.eliminate([opd,x,y,z,l,m,n,ux,uy,uz],delt,thresh = 1e-8,eliminate = eliminate)
