import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
from prtp.FlatComponent import FlatComponent
import prtp.transformationsf as trans

class CollimatorPlate(FlatComponent):
    
    ## Initialization Functions:
    
    def __init__(self,x=0,y=0,z=0,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0, l=None,w=None,collfunc=None):
        '''
        Initializes a CollimatorPlate Object, requires the following arguments:
        
        x,y,z - The position of a point along the plate
        nx,ny,nz - The components of a vector normal to the plate's surface
        sx,sy,sz - The surface vector of this plate, should be defined in a direction that makes sense for the system (straight up, for example)
        l - The length of the plate, that is, its extent in the s direction
        w - The width of the plate, that is, its extent in the sxn direction
        collfunc - A user-defined function to determine if Rays have missed the Component
        '''
        FlatComponent.__init__(self,x,y,z,nx,ny,nz,sx,sy,sz, collfunc=collfunc)
        self.l = l
        self.w = w
    
    
    def copy(self):
        '''
        Function copy:
        Returns a copy of this CollimatorPlate
        
        Inputs:
        None
        
        Outputs:
        An identical CollimatorPlate Object with the same attributes
        '''
        return CollimatorPlate(self.x,self.y,self.z,
        self.nx,self.ny,self.nz,
        self.fx,self.fy,self.fz,
        self.l,self.w,self.collisionfunction)
    
    
    ## Sample Collision Functions:
    
    def wires(self,rays,thickness=.1,sep=1):
        '''
        Function wires:
        A collision function, meant to be assigned to self.collfunc
        Simulates wires of a given thickness running infinitely in the self.surface() direction.
        The wires are separated by a distance sep
        If a photon collides with a wire, it is removed
        
        Inputs:
        rays - The rays, traced to the plate's surface, that you want to analyze
        thickness - The thickness of the wires
        sep - The distance between the centers of adjacent wires
        
        Outputs:
        A trutharray, a photon's value is True if it collides with a wire
        
        Notes:
        - Since the wires run in the self.surface() direction, only the x-position of the photons is considered
        - The "middle" wire is centered at x=0. Thus, a photon with x=0 will always be removed.
        '''
        x,y = self.getPosns(rays)
        
        # Condense the x-positions so they all exist on 0 < x < sep
        x = np.abs(x) % sep
        
        # Return True if the photon hits the wire to its left or right
        return np.logical_or((x < thickness/2),(x > sep - (thickness/2)))
        
    
    ## Removing Rays
    # Removing Rays is what a collimator does best
    
    def hit(self, rays, **kwargs):
        '''
        Function hit:
        This function simply return those rays that would be removed by this collimator plate without actually removing them
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        
        Outputs:
        tarray - A trutharray, containing True if the photon has hit the plate (and should be removed), contains False if the photon passes through the plate (and should not be removed)
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        - The function can use both a rectangle (from self.l and self.w) and a collision function
        '''
        collfuncarray = np.zeros(len(rays))
        rectarray = np.zeros(len(rays))
        # Call the collision function (if it exists)
        if self.collfunc is not None:
            collfuncarray = self.collfunc(rays,**kwargs)
        # Check if rays miss the plate (if length and width are defined)
        if (self.l is not None) and (self.w is not None):
            x,y = self.getPosns(rays)
            rectarray = np.logical_or(np.abs(x) > self.w/2, np.abs(y) > self.l/2)
        # Return an array showing if photons should be removed by collfunc or missed the plate altogether
        return np.logical_or(collfuncarray,rectarray)
    
    def removemissed(self,rays, **kwargs):
        '''
        Function removemissed:
        Removes the rays which have missed the collimator.
        This function is slightly different than Grating.removemissed.
        The Grating version will handle either the Grating's collision function or the length and width, but not both.
        CollimatorPlate.removemissed can handle both a collision function and a length/width.
        This is done to make the collision function simpler.
        Collimators will already have complicated collision functions, so having removemissed automatically handle length and width makes it much easier.
        
        Inputs:
        rays - a Rays Object which has been traced to this CollimatorPlate
        **kwargs - Any additional arguments that need to be passed to self.collfunc
        
        Output:
        Up to two tuples containing how many rays were passed in and how many survived.
        One tuple will show how many actually fell onto the rectangle of the plate (if length and width are not None)
        The other tuple will show how many photons survived the collisionfunction (if it is not None)
        
        Notes:
        -self.collfunc should return True if the photon collides with the collimator plate and needs to be removed
        '''
        # Store the length of the incoming rays
        l = len(rays)
        
        # Initialize the tuples we will eventually return
        # t1 holds the info about rays that missed the collimator plate
        # t2 holds info about rays that failed the collision function
        t1 = None
        t2 = None
        
        # get the positions of the rays:
        x,y = self.getPosns(rays)
        
        # Check if the photons missed the plate (if length and width are defined)
        if (self.l is not None) and (self.w is not None):
            tarray = np.logical_or(np.abs(x) > self.w/2, np.abs(y) > self.l/2)
            rays.remove(tarray)
            
            # Store the tuple we will eventually return
            t1 = (l,len(rays))
        
        # Get the new length of the rays
        l = len(rays)
        
        # apply the collision function if it exists
        if (self.collfunc is not None):
            tarray = self.collfunc(rays,**kwargs)
            
            rays.remove(tarray)
            
            # Store the tuple we will eventually return
            t2 = (l,len(rays))
        
        return t1,t2