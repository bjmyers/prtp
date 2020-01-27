import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
import prtp.transformationsf as trans
import astropy.units as u

class FlatComponent:
    
    ## Initialization Function
    
    @u.quantity_input(x=u.mm,y=u.mm,z=u.mm)
    def __init__(self,x=0*u.mm,y=0*u.mm,z=0*u.mm,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0, collfunc=None):
        '''
        Initializes a FlatComponent Object, requires the following arguments:
        
        x,y,z - The position of a point along the component, must be astropy
            units of length
        nx,ny,nz - The components of a vector normal to the component surface
        fx,fy,fz - The components of a direction vector along the component's 
            surface
        collfunc - A user-defined function to determine if Rays have missed the 
            Component
        
        Notes:
        - Raises an error if the normal and surface vectors are not orthogonal
        - Surface and Normal vectors can have magnitudes other than one, but 
            this function will automatically normalize them
        '''
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
        
        self.sx = sx
        self.sy = sy
        self.sz = sz
        surflen = np.sqrt(sx**2 + sy**2 + sz**2)
        self.sx /= surflen
        self.sy /= surflen
        self.sz /= surflen
        
        self.collisionfunction = collfunc
        
        # Vectors can be orthogonal within some tolerance (here 1e-8)
        dot = np.abs(np.dot([self.sx,self.sy,self.sz],[self.nx,self.ny,self.nz]))
        if np.abs(np.dot([self.sx,self.sy,self.sz],[self.nx,self.ny,self.nz])) > 1e-8:
            raise ValueError('Normal and Surface Vectors are not Orthogonal, dot product is ' + str(dot))
    
    
    ## Access Functions:
    # Used to get vectors in easy to use forms
    
    def Normal(self):
        '''
        Function Normal:
        Returns the normal Vector in a 3-element Numpy Array
        
        Inputs:
        None
        
        Outputs:
        - The normal vector of this Component in x, y, and z parts
        '''
        return np.array([self.nx,self.ny,self.nz])
    
    def Surface(self):
        '''
        Function Surface:
        Returns the surface Vector in a 3-element Numpy Array
        
        Inputs:
        None
        
        Outputs:
        - The surface vector of this Component in x, y, and z parts
        '''
        return np.array([self.sx,self.sy,self.sz])


    
    ## Spatial Manipulation Functions:
    # Functions that move and rotate the component in space
    
    @u.quantity_input(dx=u.mm,dy=u.mm,dz=u.mm)
    def translate(self,dx=0*u.mm,dy=0*u.mm,dz=0*u.mm):
        '''
        Function translate:
        Moves the Component in x,y and z
        
        Inputs:
        dx,dy,dz - The amount to move in x, y, and z, respectively, must be 
            astropy units of length
        
        Outputs:
        None
        
        Notes:
        - This move is relative, not absolute. That is, you will move BY dx, dy, 
            and z, you will not move TO dx, dy, and dz
        '''
        self.x += dx.to(u.mm)
        self.y += dy.to(u.mm)
        self.z += dz.to(u.mm)
    
    @u.quantity_input(theta=u.rad)
    def unitrotate(self,theta=0*u.rad,axis=1):
        '''
        Function unitrotate:
        Rotates the Component about one of the three unit vectors (x, y, or z) 
            by an amount theta
        
        Inputs:
        theta - The angle by which you want to rotate, must be an astropy unit
            of angle
        axis - integer input of 1, 2, or 3 to rotate about the x, y, or z axes, 
            respectively.
        
        Outputs:
        None
        '''
        theta = theta.to(u.rad).value
        self.nx,self.ny,self.nz = trans.rotatevector(self.nx,self.ny,self.nz,theta,axis)
        self.sx,self.sy,self.sz = trans.rotatevector(self.sx,self.sy,self.sz,theta,axis)
    
    @u.quantity_input(theta=u.rad)
    def rotate(self,theta=0*u.rad,ux=1,uy=0,uz=0):
        '''
        Function rotate:
        Rotates the Component about an arbitrary axis, defined by ux, uy, and uz
        
        Inputs:
        theta - The angle by which you want to rotate, must be an
            astropy unit of angle
        ux,uy,uz - The x, y, and z components of the vector about which you want
            to rotate the Component
        
        Outputs:
        None
        '''
        theta = theta.to(u.rad).value
        self.nx,self.ny,self.nz,tx,ty,tz = trans.rotateaxis(self.nx,self.ny,self.nz,ux,uy,uz,theta)
        self.sx,self.sy,self.sz,tx,ty,tz = trans.rotateaxis(self.sx,self.sy,self.sz,ux,uy,uz,theta)
    
    @u.quantity_input(theta=u.rad)
    def pitch(self,theta=0*u.rad):
        '''
        Function pitch:
        Rotates the Component about its sxn vector in a right-handed fashion
        
        Inputs:
        theta - The angle by which you want to rotate, must be an astropy unit
            of angle
        
        Outputs:
        None
        '''
        theta = theta.to(u.rad).value
        sxn = np.cross(self.Surface(),self.Normal())
        
        self.nx,self.ny,self.nz,tx,ty,tz = trans.rotateaxis(self.nx,self.ny,self.nz,sxn[0],sxn[1],sxn[2],theta)
        self.sx,self.sy,self.sz,tx,ty,tz = trans.rotateaxis(self.sx,self.sy,self.sz,sxn[0],sxn[1],sxn[2],theta)
    
    @u.quantity_input(theta=u.rad)
    def roll(self,theta=0*u.rad):
        '''
        Function pitch:
        Rotates the Component about its surface vector in a right-handed fashion
        
        Inputs:
        theta - The angle by which you want to rotate, must be an astropy unit
            of angle
        
        Outputs:
        None
        '''
        theta = theta.to(u.rad).value
        surf = self.Surface()
        
        self.nx,self.ny,self.nz,tx,ty,tz = trans.rotateaxis(self.nx,self.ny,self.nz,surf[0],surf[1],surf[2],theta)
        self.sx,self.sy,self.sz,tx,ty,tz = trans.rotateaxis(self.sx,self.sy,self.sz,surf[0],surf[1],surf[2],theta)
    
    @u.quantity_input(theta=u.rad)
    def yaw(self,theta=0*u.rad):
        '''
        Function pitch:
        Rotates the Component about its normal vector in a right-handed fashion
        
        Inputs:
        theta - The angle by which you want to rotate, must be an astropy unit
            of angle
        
        Outputs:
        None
        '''
        theta = theta.to(u.rad).value
        norm = self.Normal()
        
        self.nx,self.ny,self.nz,tx,ty,tz = trans.rotateaxis(self.nx,self.ny,self.nz,norm[0],norm[1],norm[2],theta)
        self.sx,self.sy,self.sz,tx,ty,tz = trans.rotateaxis(self.sx,self.sy,self.sz,norm[0],norm[1],norm[2],theta)



    
    ## Ray Tracing Functions:
    
    def trace_to_surf(self, rays, modify=True):
        '''
        Function trace_to_surf:
        Given a Rays object, traces the rays to the surface of the component
        
        Inputs:
        rays - The Rays object you want to trace
        modify - Boolean. If True, the original Rays object will be modified 
            into one on the Component's surface. If False, a copy of the Rays 
            object will be made and returned in its modified form, leaving the 
            original unchanged.
        
        Outputs:
        Rays - The Rays object with photons traced to the Component's Surface
        
        Notes:
        - Rays object will have NaNs in its x,y, and z positions if it is 
            parallel to the plane
        - A Rays object will always be returned, but if modify=False, it will 
            not be the same Rays object as the one input
        '''
        if (not modify):
            rays = rays.copy()
        
        # Calculate the distance each ray needs to travel to reach the plane
        # This distance is {deltap (dot) normal} / {direction (dot) normal}
        # where deltap is the vector pointing from the ray to the center
        
        dist = (((self.x.value - rays.x) * self.nx) + ((self.y.value - rays.y) * self.ny) + ((self.z.value - rays.z) * self.nz)) / ((rays.l * self.nx) + (rays.m * self.ny) + (rays.n * self.nz))
        
        # Move the rays that distance
        vel = np.sqrt(rays.l**2 + rays.m**2 + rays.n**2)
        rays.set(x = rays.x + dist*rays.l/vel,y = rays.y + dist*rays.m/vel,z = rays.z + dist*rays.n/vel)
        
        # Update surface normal vectors
        rays.move(ux = self.nx, uy = self.ny, uz = self.nz)
        
        return rays
    
    def getPosns(self,rays):
        '''
        Function getPosns:
        Given Rays that have been traced to the Component, finds their x and y 
            positions on the Plane. The X position is its distance from the 
            center in the sxn direction. The Y position is its distance from the 
            center in the s direction. Therefore, the n axis is the z-direction 
            basis vector.
        
        Inputs:
        rays - A Rays object that has been traced to the Component Plane
        
        Outputs:
        x,y - Arrays holding the X and Y position of each photon
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up 
            to the user. It will still create output, but it will be wonky
        '''
        if (len(rays) == 0):
            return np.array([]),np.array([])
        
        
        # Define the x-direction basis vector, the cross product of s and n
        xdir = np.cross((self.sx,self.sy,self.sz),(self.nx,self.ny,self.nz))
        ydir = np.array([self.sx,self.sy,self.sz])
        
        # X-components are {deltap (dot) (s x n)}
        xs = (rays.x - self.x.value) * (xdir[0]) + (rays.y - self.y.value) * (xdir[1]) + (rays.z - self.z.value) * (xdir[2])
        
        # Y-components are {deltap (dot) s}
        ys = (rays.x - self.x.value) * (ydir[0]) + (rays.y - self.y.value) * (ydir[1]) + (rays.z - self.z.value) * (ydir[2])
        
        return xs,ys
    
    def getDist(self,rays):
        '''
        Function getDist:
        Given a Rays object, this function finds the distance from each Ray to
            the plane of this component
        
        Inputs:
        rays - A rays object
        
        Outputs:
        An array of length len(rays) containing the distance of each photon to
            the plane
        '''
        # Calculate the distance each ray needs to travel to reach the plane
        # This distance is {deltap (dot) normal} / {direction (dot) normal}
        # where deltap is the vector pointing from the ray to the center
        
        return (((self.x.value - rays.x) * self.nx) + ((self.y.value - rays.y) * self.ny) + ((self.z.value - rays.z) * self.nz)) / ((rays.l * self.nx) + (rays.m * self.ny) + (rays.n * self.nz))
    
    
    def reflect(self,rays):
        '''
        Function reflect:
        Given Rays that have been traced to the Component, reflectes the Rays 
            off of the surface.
        
        Inputs:
        rays - A Rays object that has been traced to the Component Plane
        
        Outputs:
        Nothing
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up 
            to the user.
        '''
        length = len(rays)
        rays.ux = np.ones(length) * self.nx
        rays.uy = np.ones(length) * self.ny
        rays.uz = np.ones(length) * self.nz
        rays.reflect()