import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
import prtp.specialFunctions as func
import prtp.transformationsf as trans
import astropy.units as u

class WolterOptic:
    
    #TODO: Rays are potentially being traced to the wrong side. This should be a focus of future testing
    
    def __init__(self,x=0*u.mm,y=0*u.mm,z=0*u.mm,nx=0,ny=0,nz=1,
                r0=1*u.mm,z0=1*u.mm,psi=1):
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
        '''
        if (type(x) != u.quantity.Quantity or type(y) != u.quantity.Quantity 
        or type(z) != u.quantity.Quantity or type(r0) != u.quantity.Quantity
        or type(z0) != u.quantity.Quantity):
            raise ValueError('x, y, z, r0, and z0 must all be astropy units of length')
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
