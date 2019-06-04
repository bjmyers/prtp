import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
import prtp.specialFunctions as func
import prtp.transformationsf as trans
from prtp.Combination import Combination
import astropy.units as u


class WolterOptic:
    
    def __init__(self,x=0*u.mm,y=0*u.mm,z=0*u.mm,nx=0,ny=0,nz=1,
                r0=1*u.mm,z0=1*u.mm,psi=1,axial_length = None,mirror_sep = None):
        '''
        WolterOptic Object:
        
        Parameters:
        x,y,z - The Cartesian Coordinates of the focus of the optic, must be 
            astropy units of length
        nx,ny,nz - The components of the vector pointing outwards from the focus 
            of the optic
        r0 - The radius of the optics when they converge, must be an astropy 
            unit of length
        z0 - The position along the n-direction at which the optics converge 
            (the focus is at position 0) Must be an astropy unit of length
        psi - Some ratio of angles, I don't know, I think its always 1
        axial_length - The length along the n-axis of the mirror, must be an 
            astropy unit of length
        mirror_sep - The separation between two mirrors, must be an astropy unit
            of length
        
        Notes:
        - axial_length and mirror_sep are defined such that the primary mirror
            will run from:
            [z0 + (mirror_sep/2)] to [z0 + (mirror_sep/2) + axial_length]
            and the secondary mirror will run from:
            [z0 - (mirror_sep/2)] to [z0 - (mirror_sep/2) - axial_length]
        - If axial_length and mirror_sep are defined, any photons which miss
            the mirror in the n-direction will removed.
        '''
        if (type(x) != u.quantity.Quantity or type(y) != u.quantity.Quantity 
        or type(z) != u.quantity.Quantity or type(r0) != u.quantity.Quantity
        or type(z0) != u.quantity.Quantity):
            raise ValueError('x, y, z, r0, and z0 must all be astropy units of length')
            
        if (axial_length is not None or mirror_sep is not None):
            if (type(axial_length) != u.quantity.Quantity or 
                type(mirror_sep) != u.quantity.Quantity):
                raise ValueError('axial_length and mirror_sep must both be defined and must be astropy units of length')
            else:
                self.axial_length = axial_length.to(u.mm)
                self.mirror_sep = mirror_sep.to(u.mm)
        else:
            self.axial_length = None
            self.mirror_sep = None
            
        self.x = x.to(u.mm)
        self.y = y.to(u.mm)
        self.z = z.to(u.mm)
        self.nx = nx
        self.ny = ny
        self.nz = nz
        
        normlen = np.sqrt(nx**2 + ny**2 + nz**2)
        self.nx /= normlen
        self.ny /= normlen
        self.nz /= normlen
        
        self.r0 = r0.to(u.mm)
        self.z0 = z0.to(u.mm)
        self.psi = psi
    
    def copy(self):
        '''
        Function copy:
        Returns a copy of this Wolter Optic
        
        Inputs:
        None
        
        Outputs:
        An identical WolterOptic Object with the same attributes
        '''
        return WolterOptic(self.x,self.y,self.z,
        self.nx,self.ny,self.nz,
        self.r0,self.z0,self.psi)
    
    
    ## Spatial Manipulation Functions:
    
    def translate(self,dx,dy,dz):
        '''
        Function translate:
        Moves the Optic in x,y and z
        
        Inputs:
        dx,dy,dz - The amount to move in x, y, and z, respectively. Must be 
            astropy units of length
        
        Outputs:
        None
        
        Notes:
        - This move is relative, not absolute. That is, you will move BY dx, dy,
            and z, you will not move TO dx, dy, and dz
        '''
        if (type(dx) != u.quantity.Quantity or type(dy) != u.quantity.Quantity
        or type(dz) != u.quantity.Quantity):
            raise ValueError('dx, dy, and dz must be astropy units of length')
        self.x += dx.to(u.mm).value
        self.y += dy.to(u.mm).value
        self.z += dz.to(u.mm).value
    
    def unitrotate(self,theta,axis):
        '''
        Function unitrotate:
        Rotates the Optic about one of the three unit vectors (x, y, or z) by an 
            amount theta
        
        Inputs:
        theta - The angle by which you want to rotate, in radians. Must be an
            astropy unit of angle
        axis - integer input of 1, 2, or 3 to rotate about the x, y, or z axes, 
            respectively.
        
        Outputs:
        None
        '''
        if (type(theta) != u.quantity.Quantity):
            raise ValueError('theta must be an astropy unit of angle')
        self.nx,self.ny,self.nz = trans.rotatevector(self.nx,self.ny,self.nz,theta.to(u.rad).value,axis)
    
    def rotate(self,theta,ux,uy,uz):
        '''
        Function rotate:
        Rotates the Optic about an arbitrary axis, defined by ux, uy, and uz
        
        Inputs:
        theta - The angle by which you want to rotate, in radians. Must be and
            astropy unit of angle
        ux,uy,uz - The x, y, and z components of the vector about which you want 
            to rotate the Optic
        
        Outputs:
        None
        '''
        if (type(theta) != u.quantity.Quantity):
            raise ValueError('theta must be an astropy unit of angle')
        self.nx,self.ny,self.nz,tx,ty,tz = trans.rotateaxis(self.nx,self.ny,self.nz,ux,uy,uz,theta.to(u.rad).value)


    ## Tracing Rays to the Optic:
    
    def tracehelper(self,rays,func,autoreflect=True,considerweights=False):
        '''
        Function tracehelper:
        Allows Wolter Optics to trace their rays to their surface
        
        Inputs:
        rays - The Rays object you want to trace to the optic
        func - A The function used to trace the photons to a surface, this is 
            handled by each subclass
        autoreflect - If True, photons will be reflected off of the optic's 
            surface automatically
        considerweights - Should be True if the photons are weighted
        
        Outputs:
        A tuple containing the original number of rays and the final number of 
            rays (that hit the optic)
        '''
        
        # First find the angles through which we need to rotate the rays
        norm = np.array([self.nx,self.ny,self.nz])
        
        # Check if optic is already aligned, if it is our angle-finding code will produce an error, so we need to skip it
        if norm[0] == 0 and norm[1] == 0:
            # If normal is antiparallel with the z-axis, we need to rotate 180 degrees about x
            if np.sign(norm[2]) < 0:
                thetax = np.pi
            else:
                thetax = 0
            thetaz = 0
        else:
        
            # Find the angle between the normal vector projected to the xy-plane and the y-axis, this is how far we will need to rotate our velocities about the z-axis
            # [0,1] is the unit vector in the y-direction
            thetaz = np.arccos(np.dot([norm[0],norm[1]],[0,1])/np.sqrt(norm[0]**2 + norm[1]**2))
            
            # Check if normal vector projected to xy-plane is in Quadrants 2 or 3, if so, we need to change the sign of our angle
            if (norm[0] < 0):
                thetaz *= -1
            
            # Rotate the normal vector so its in the yz-plane
            norm = trans.rotatevector(norm[0],norm[1],norm[2],thetaz,3)
            
            # Now we need to do the same thing in the yz-plane. Find the angle we need to rotate about x in order to align our normal vector with the z-axis
            # In this plane, [0,1] is the unit vector in the z-direction
            thetax = np.arccos(np.dot([norm[1],norm[2]],[0,1])/np.sqrt(norm[1]**2 + norm[2]**2))
            
            # Check if normal vector projected to xy-plane is in Quadrants 2 or 3, if so, we need to change the sign of our angle
            if (norm[1] < 0):
                thetax *= -1
        
        # Translates the rays so that the focus of the optic is at the origin
        rays.translate(-self.x.value,-self.y.value,-self.z.value)
        
        # Rotates the rays so they're aligned with the optic
        rays.rotatevector('pos',thetaz,3)
        rays.rotatevector('pos',thetax,1)
        rays.rotatevector('dir',thetaz,3)
        rays.rotatevector('dir',thetax,1)
        rays.rotatevector('norm',thetaz,3)
        rays.rotatevector('norm',thetax,1)
        
        # At this point, the photons have been transformed, but each type of Wolter Optic has a different way of tracing rays.
        func(rays,autoreflect)
        
        # Undo all of our previous transformations
        rays.rotatevector('pos',-thetax,1)
        rays.rotatevector('pos',-thetaz,3)
        rays.rotatevector('dir',-thetax,1)
        rays.rotatevector('dir',-thetaz,3)
        rays.rotatevector('norm',-thetax,1)
        rays.rotatevector('norm',-thetaz,3)
        rays.translate(self.x.value,self.y.value,self.z.value)
        
        # Find how many rays missed the optic
        inputlength = rays.length(considerweights)
        
        # Remove missed photons
        rays.remove(np.isnan(rays.x))
        
        # Return the efficiency of the optic
        return ("Missed Optic",inputlength,rays.length(considerweights))


class WolterPrimary(WolterOptic):
    
    def __init__(self,x=0*u.mm,y=0*u.mm,z=0*u.mm,nx=0,ny=0,nz=1,
        r0=1*u.mm,z0=1*u.mm,psi=1,axial_length = None,mirror_sep = None):
        '''
        WolterPrimary Object:
        
        Parameters:
        x,y,z - The Cartesian Coordinates of the focus of the optic. Must be 
            astropy units of length
        nx,ny,nz - The components of the vector pointing outwards from the focus 
            of the optic
        r0 - The radius of the optic when it converges with the WolterSecondary.
            Must be an astropy unit of length
        z0 - The position along the n-direction at which this optic converges,
            must be an astropy unit of length
            with the WolterSecondary (the focus is at position 0)
        psi - Some ratio of angles, I don't know, I think its always 1
        axial_length - The length along the n-axis of the mirror, must be an 
            astropy unit of length
        mirror_sep - The separation between two mirrors, must be an astropy unit
            of length
        '''
        WolterOptic.__init__(self,x,y,z,nx,ny,nz,r0,z0,psi,axial_length,mirror_sep)


    ## Tracing Rays to the Optic:
    
    def tracefunction(self,rays,autoreflect=True):
        '''
        Function tracefunction:
        Determines what happens to the rays once they has been transformed, 
            should only be called by the WolterOptic superclass's trace function
        '''
        rays.wolterprimary(self.r0.value,self.z0.value,self.psi)
        
        # Consider which photons missed the primary mirror
        if self.axial_length is not None and self.mirror_sep is not None:
            tarray = np.logical_or((rays.z < (self.z0 + self.mirror_sep/2).value),
                    rays.z > (self.z0 + self.mirror_sep/2 + self.axial_length).value)
            rays.remove(tarray)
        
        if autoreflect:
            rays.reflect()
    
    def trace(self,rays,autoreflect=True,remove=True,considerweights=False):
        '''
        Function trace:
        Traces rays to the wolter optic
        
        Inputs:
        rays - The Rays object you want to trace to the optic
        autoreflect - If True, photons will be reflected off of the optic's 
            surface automatically
        considerweights - Should be True if the photons are weighted
        
        Outputs:
        eff - A tuple containing information about how many photons hit the
            optic
        '''
        l = rays.length(considerweights)
        self.tracehelper(rays,self.tracefunction,autoreflect,considerweights)
        return ("Missed Primary Optic",l,rays.length(considerweights))


class WolterSecondary(WolterOptic):
    
    def __init__(self,x=0*u.mm,y=0*u.mm,z=0*u.mm,nx=0,ny=0,nz=1,
        r0=1*u.mm,z0=1*u.mm,psi=1,axial_length = None,mirror_sep = None):
        '''
        WolterSecondary Object:
        
        Parameters:
        x,y,z - The Cartesian Coordinates of the focus of the optic. Must be
            astropy units of length
        nx,ny,nz - The components of the vector pointing outwards from the focus 
            of the optic
        r0 - The radius of the optic when it converges with the WolterPrimary.
            Must be an astropy unit of length
        z0 - The position along the n-direction at which this optic converges 
            with the WolterPrimary (the focus is at position 0). Must be an 
            astropy unit of length
        psi - Some ratio of angles, I don't know, I think its always 1
        axial_length - The length along the n-axis of the mirror, must be an 
            astropy unit of length
        mirror_sep - The separation between two mirrors, must be an astropy unit
            of length
        '''
        WolterOptic.__init__(self,x,y,z,nx,ny,nz,r0,z0,psi,axial_length,mirror_sep)


    ## Tracing Rays to the Optic:
    
    def tracefunction(self,rays,autoreflect=True):
        '''
        Function tracefunction:
        Determines what happens to the rays once they has been transformed, 
            should only be called by the WolterOptic superclass's trace function
        '''
        rays.woltersecondary(self.r0.value,self.z0.value,self.psi)
        
        # Consider which photons missed the mirror
        if self.axial_length is not None and self.mirror_sep is not None:
            tarray = np.logical_or((rays.z > (self.z0 - self.mirror_sep/2).value),
                    rays.z < (self.z0 - self.mirror_sep/2 - self.axial_length).value)
            rays.remove(tarray)
        
        if autoreflect:
            rays.reflect()
    
    def trace(self,rays,autoreflect=True,considerweights=False):
        '''
        Function trace:
        Traces rays to the wolter optic
        
        Inputs:
        rays - The Rays object you want to trace to the optic
        autoreflect - If True, photons will be reflected off of the optic's 
            surface automatically
        considerweights - Should be True if the photons are weighted
        
        Outputs:
        eff - A tuple containing information about how many photons hit the
            optic
        '''
        l = rays.length(considerweights)
        self.tracehelper(rays,self.tracefunction,autoreflect,considerweights)
        return ("Missed Secondary Optic",l,rays.length(considerweights))


class WolterTypeOne(WolterOptic):
    
    def __init__(self,x=0*u.mm,y=0*u.mm,z=0*u.mm,nx=0,ny=0,nz=1,
    r0=1*u.mm,z0=1*u.mm,psi=1,axial_length = None,mirror_sep = None):
        '''
        WolterTypeOne Object:
        - Combines a WolterPrimary and a WolterSecondary Object into one single 
            object
        The WolterTypeOne Object will process rays faster than you would with 
            the two individual objects, owing to the fact that the rays only 
                have to be transformed once.
        DO NOT use a WolterTypeOne object if:
            - You need to move the two components relative to each other so that 
                they no longer share a common focus (e.g: you try to handle 
                misalignments in the two optics separately)
            - You need to perform a tracing more complicated than:
                - Trace to Primary
                - Reflect
                - Trace to Secondary
            e.g: with WolterTypeOne Objects you cannot add scattering on the 
                surface of the Primary (since the rays are automatically traced 
                to the secondary, giving the user no time to perform further 
                modifications)
        
        
        Parameters:
        x,y,z - The Cartesian Coordinates of the focus of the optic. Must be
            astropy units of length
        nx,ny,nz - The components of the vector pointing outwards from the focus 
            of the optic
        r0 - The radius of the optic when it converges with the WolterSecondary.
            Must be an astropy unit of length
        z0 - The position along the n-direction at which this optic converges .
            Must be an astropy unit of length
            with the WolterSecondary (the focus is at position 0)
        psi - Some ratio of angles, I don't know, I think its always 1
        axial_length - The length along the n-axis of the mirror, must be an 
            astropy unit of length
        mirror_sep - The separation between two mirrors, must be an astropy unit
            of length
        '''
        WolterOptic.__init__(self,x,y,z,nx,ny,nz,r0,z0,psi,axial_length,mirror_sep)


    ## Tracing Rays to the Optic:
    
    def tracefunction(self,rays,autoreflect,eliminate='remove'):
        '''
        Function tracefunction:
        Determines what happens to the rays once they has been transformed, 
            should only be called by the WolterOptic superclass's trace function
        
        Notes:
        -The eliminate parameter has two options, "remove" and "nan". For almost
            all cases it will be "remove." But WolterModule objects require 
            it to be set to "nan", in which eliminated rays are set to NaN,
            rather than removed altogether from the Rays object
        '''
        rays.wolterprimary(self.r0.value,self.z0.value,self.psi)
        
        # Consider which photons missed the primary mirror
        if self.axial_length is not None and self.mirror_sep is not None:
            tarray = np.logical_or((rays.z < (self.z0 + self.mirror_sep/2).value),
                    rays.z > (self.z0 + self.mirror_sep/2 + self.axial_length).value)
            if (eliminate.lower() == 'remove'):
                rays.remove(tarray)
            else:
                rays.x[tarray] = np.nan
        
        rays.reflect()
        rays.woltersecondary(self.r0.value,self.z0.value,self.psi)
        
        # When called with eliminate='nan', this block will cause runtime errors,
        # we will ignore them using this next line
        with np.errstate(invalid='ignore'):
        
            # Consider which photons missed the secondary mirror
            if self.axial_length is not None and self.mirror_sep is not None:
                tarray = np.logical_or((rays.z > (self.z0 - self.mirror_sep/2).value),
                        rays.z < (self.z0 - self.mirror_sep/2 - self.axial_length).value)
                if (eliminate.lower() == 'remove'):
                    rays.remove(tarray)
                else:
                    rays.x[tarray] = np.nan
        
        if autoreflect:
            rays.reflect()
    
    def trace(self,rays,autoreflect=True,considerweights=False):
        '''
        Function trace:
        Traces rays to the wolter optic
        
        Inputs:
        rays - The Rays object you want to trace to the optic
        autoreflect - If True, photons will be reflected off of the optic's 
            surface automatically
        considerweights - Should be True if the photons are weighted
        
        Outputs:
        eff - A tuple containing information about how many photons hit the
            optic
        '''
        l = rays.length(considerweights)
        self.tracehelper(rays,self.tracefunction,autoreflect,considerweights)
        return ("Missed Wolter Optic",l,rays.length(considerweights))


class WolterModule(WolterOptic, Combination):
    
    @u.quantity_input(x=u.mm,y=u.mm,z=u.mm,r0=u.mm,z0=u.mm,axial_length=u.mm,mirror_sep=u.mm)
    def __init__(self, x=0*u.mm,y=0*u.mm,z=0*u.mm,nx=0,ny=0,nz=1, 
                r0=None, z0 = None, psi = 1,axial_length = None,mirror_sep = None):
        '''
        Function __init__:
        Initializes a WolterModule Object. This object consists of a list of 
            WolterTypeOne objects in the same way that a GratingStack consists
            of a list of Grating Objects. The objects in a WolterModule Object
            are nested together.
        
        Inputs:
        x,y,z - The position of the focus of the optic, must be astropy units
            of length
        nx,ny,nz - Components of a normal vector pointing outwards from the focus 
            of the optic
        r0 - The radius of the optic when it converges with the WolterSecondary.
            Must be an astropy unit of length
        z0 - The position along the n-direction at which this optic converges .
            Must be an astropy unit of length
            with the WolterSecondary (the focus is at position 0)
        psi - Some ratio of angles, I don't know, I think its always 1
        axial_length - The length along the n-axis of the mirror, must be an 
            astropy unit of length
        mirror_sep - The separation between two mirrors, must be an astropy unit
            of length
        
        
        Notes:
        - r0,z0,axial_length,and mirror_sep can be lists or numpy arrays of the 
            same length with astropy units of length attached. If they are lists
            of length n, then n WolterTypeOne objects will be added to the 
            componentlist of this Module. Otherwise, Mirrors have to be added
            individually.
        - WolterModule also descends from WolterOptic, and as such has its own
            parameters for x,y,z,nx,ny, and nz. Since every mirror in the Module
            must have the same value for x,y,z,nx,ny, and nz, the Module's values
            for these parameters is just the same as the values for each mirror.
            These parameters are needed for the orientation algorithm 
            tracehelper()
        '''
        Combination.__init__(self)
        WolterOptic.__init__(self,x,y,z,nx,ny,nz)
        if (r0 is not None or z0 is not None):
            try:
                lr = len(r0)
                lz = len(z0)
                la = len(axial_length)
                lm = len(mirror_sep)
                if (lr != lz or lr != la or lr != lm):
                    raise ValueError('r0, z0, axial_length, and mirror_sep must be the same length!')
                # Iterate through each list and add the Wolter Optic to 
                # the component list
                for i in range(lr):
                    self.componentlist.append(WolterTypeOne(x,y,z,nx,ny,nz,r0[i],z0[i],psi,axial_length[i],mirror_sep[i]))
            except:
                # This will execute if r0 and z0 are single values
                self.componentlist.append(WolterTypeOne(x,y,z,nx,ny,nz,r0,z0,psi,axial_length,mirror_sep))
    
    ## Tracing Rays to the Optic:
    def tracefunction(self,rays,autoreflect=True):
        '''
        Function defaultTrace:
        Traces the Rays through the Mirror Module in the order that the Mirrors
            were added to the Module.
        
        Inputs:
        rays - The rays you want to trace through the stack
        autoreflect - Boolean saying if you want to automatically reflect the
            rays.
        
        Outputs:
        Nothing
        '''
        # Make a blank Rays object to store the Rays that make it
        finalrays = Rays()
        
        # Keep track of the input rays for when we're finished with one Mirror
        inputrays = rays.copy()
        
        # Iterate through each Mirror Object
        for r in self.componentlist:
            # Through each pass we need to ensure that the rays that make it are 
            # placed into a final rays object
            # All those that miss are passed to the next Mirror
            r.tracefunction(rays,autoreflect,eliminate='nan')

            # Find the Rays which missed the 
            tarray = np.logical_not(np.isnan(rays.x))
            hitrays = rays.split(tarray)
            
            # Make sure at least some rays have hit the grating
            if (len(hitrays) == 0):
                rays = inputrays.copy()
                continue
            
            # Add the hitrays to our final tally
            finalrays += hitrays
            
            # Take the rays that hit this grating out of the original Rays object
            inputrays.remove(tarray)
            
            # Back remaining rays up to their original position
            rays = inputrays.copy()
            
            # If there are no rays left, we can stop
            if len(rays) == 0:
                break
        
        # Make it so that the original rays now contain the output
        rays.makecopy(finalrays)
    
    
    def trace(self,rays,autoreflect=True,considerweights=False):
        '''
        Function trace:
        Traces rays to the wolter module
        
        Inputs:
        rays - The Rays object you want to trace to the optic
        autoreflect - If True, photons will be reflected off of the optic's 
            surface automatically
        considerweights - Should be True if the photons are weighted
        
        Outputs:
        eff - A tuple containing information about how many photons hit the
            optic
        '''
        l = rays.length(considerweights)
        self.tracehelper(rays,self.tracefunction,autoreflect,considerweights)
        return ("Missed Wolter Module",l,rays.length(considerweights))