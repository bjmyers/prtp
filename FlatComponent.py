import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
import prtp.transformationsf as trans

class FlatComponent:
    
    ## Initialization Function
    
    def __init__(self,x=0,y=0,z=0,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0, collfunc=None):
        '''
        Initializes a FlatComponent Object, requires the following arguments:
        
        x,y,z - The position of a point along the component
        nx,ny,nz - The components of a vector normal to the component surface
        fx,fy,fz - The components of a direction vector along the component's surface
        collfunc - A user-defined function to determine if Rays have missed the Component
        
        Notes:
        - Raises an error if the normal and surface vectors are not orthogonal
        - Surface and Normal vectors can have magnitudes other than one, but this function will automatically normalize them
        '''
        self.x = x
        self.y = y
        self.z = z
        
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
        if np.abs(np.dot([self.sx,self.sy,self.sz],[self.nx,self.ny,self.nz])) > 1e-8:
            raise ValueError('Normal and Surface Vectors are not Orthogonal')
    
    
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
    
    def translate(self,dx,dy,dz):
        '''
        Function translate:
        Moves the Component in x,y and z
        
        Inputs:
        dx,dy,dz - The amount to move in x, y, and z, respectively
        
        Outputs:
        None
        
        Notes:
        - This move is relative, not absolute. That is, you will move BY dx, dy, and z, you will not move TO dx, dy, and dz
        '''
        self.x += dx
        self.y += dy
        self.z += dz
    
    def unitrotate(self,theta,axis):
        '''
        Function unitrotate:
        Rotates the Component about one of the three unit vectors (x, y, or z) by an amount theta
        
        Inputs:
        theta - The angle by which you want to rotate, in radians
        axis - integer input of 1, 2, or 3 to rotate about the x, y, or z axes, respectively.
        
        Outputs:
        None
        '''
        self.nx,self.ny,self.nz = trans.rotatevector(self.nx,self.ny,self.nz,theta,axis)
        self.sx,self.sy,self.sz = trans.rotatevector(self.sx,self.sy,self.sz,theta,axis)
    
    def rotate(self,theta,ux,uy,uz):
        '''
        Function rotate:
        Rotates the Component about an arbitrary axis, defined by ux, uy, and uz
        
        Inputs:
        theta - The angle by which you want to rotate, in radians
        ux,uy,uz - The x, y, and z components of the vector about which you want to rotate the Component
        
        Outputs:
        None
        '''
        self.nx,self.ny,self.nz,tx,ty,tz = trans.rotateaxis(self.nx,self.ny,self.nz,ux,uy,uz,theta)
        self.sx,self.sy,self.sz,tx,ty,tz = trans.rotateaxis(self.sx,self.sy,self.sz,ux,uy,uz,theta)



    
    ## Ray Tracing Functions:
    
    def trace_to_surf(self, rays, modify=True):
        '''
        Function trace_to_surf:
        Given a Rays object, traces the rays to the surface of the component
        
        Inputs:
        rays - The Rays object you want to trace
        modify - Boolean. If True, the original Rays object will be modified into one on the Component's surface. If False, a copy of the Rays object will be made and returned in its modified form, leaving the original unchanged.
        
        Outputs:
        Rays - The Rays object with photons traced to the Component's Surface
        
        Notes:
        - Rays object will have NaNs in its x,y, and z positions if it is parallel to the plane
        - A Rays object will always be returned, but if modify=False, it will not be the same Rays object as the one input
        '''
        if (not modify):
            rays = rays.copy()
        
        gratx = np.ones(len(rays))*self.x
        graty = np.ones(len(rays))*self.y
        gratz = np.ones(len(rays))*self.z
        gratnx = np.ones(len(rays))*self.nx
        gratny = np.ones(len(rays))*self.ny
        gratnz = np.ones(len(rays))*self.nz
        
        diff = np.array([gratx - rays.x,graty - rays.y, gratz - rays.z])
        # Calculate the distance each ray needs to travel to reach the plane
        #Note: (x*y).sum(0) finds the dot product of two arrays of vectors
        dist = (diff*np.array([gratnx,gratny,gratnz])).sum(0) / (np.array([rays.l,rays.m,rays.n])*np.array([gratnx,gratny,gratnz])).sum(0)
        
        # Move the rays that distance
        vel = np.sqrt(rays.l**2 + rays.m**2 + rays.n**2)
        rays.set(x = rays.x + dist*rays.l/vel,y = rays.y + dist*rays.m/vel,z = rays.z + dist*rays.n/vel)
        
        return rays
    
    def getPosns(self,rays):
        '''
        Function getPosns:
        Given Rays that have been traced to the Component, finds their x and y positions on the Plane. The X position is its distance from the center in the sxn direction. The Y position is its distance from the center in the s direction. Therefore, the n axis is the z-direction basis vector.
        
        Inputs:
        rays - A Rays object that has been traced to the Component Plane
        
        Outputs:
        x,y - Arrays holding the X and Y position of each photon
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user. It will still create output, but it will be wonky
        '''
        # Define the x-direction basis vector, the cross product of s and n
        xdir = np.cross((self.sx,self.sy,self.sz),(self.nx,self.ny,self.nz))
        ydir = np.array([self.sx,self.sy,self.sz])
        
        # Stack direction vectors so they can be easily compared with the rays
        xdir = np.vstack([xdir] * len(rays))
        ydir = np.vstack([ydir] * len(rays))
        
        # Define the vector going from the rays to the center of the plane
        deltap = np.array([self.x - rays.x, self.y - rays.y, self.z - rays.z]).transpose()
        
        # X-components are {deltap (dot) (s x n)}
        xs = (deltap * xdir).sum(1)
        # Y-components are {deltap (dot) s}
        ys = (deltap * ydir).sum(1)
        
        return xs,ys

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    