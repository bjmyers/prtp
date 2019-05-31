import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
from prtp.FlatComponent import FlatComponent
import prtp.transformationsf as trans
import astropy.units as u

class CollimatorPlate(FlatComponent):
    
    ## Initialization Functions:
    
    def __init__(self,x=0,y=0,z=0,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0, l=None,w=None,collfunc=None):
        '''
        Initializes a CollimatorPlate Object, requires the following arguments:
        
        x,y,z - The position of a point along the plate. Must be astropy units
            of length
        nx,ny,nz - The components of a vector normal to the plate's surface
        sx,sy,sz - The surface vector of this plate, should be defined in a 
            direction that makes sense for the system (straight up, for example)
        l - The length of the plate, that is, its extent in the s direction, 
            must be an astropy unit of length
        w - The width of the plate, that is, its extent in the sxn direction,
            must be an astropy unit of length
        collfunc - A user-defined function to determine if Rays have missed the 
            Component
        
        Notes:
        - collfunc must take in only a Rays object. Any other argument needed
            must be given the the CollimatorPlate Object as a parameter, see
            self.wires() for an example. This restriction is done because
            an Instrument Object cannot pass any necessary arguments to 
            collfunc once it is running
        '''
        FlatComponent.__init__(self,x,y,z,nx,ny,nz,sx,sy,sz, collfunc=collfunc)
        if l is not None or w is not None:
            if (type(l) != u.quantity.Quantity or type(w) != u.quantity.Quantity):
                raise ValueError('l and w must be astropy units of length')
        self.l = l.to(u.mm)
        self.w = w.to(u.mm)
    
    
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
        self.sx,self.sy,self.sz,
        self.l,self.w,self.collisionfunction)
    
    
    ## Sample Collision Functions:
    
    def wires(self,rays):
        '''
        Function wires:
        A collision function, meant to be assigned to self.collfunc
        Simulates wires of a given thickness running infinitely in the 
            self.surface() direction.
        The wires are separated by a distance sep
        If a photon collides with a wire, it is removed
        
        Inputs:
        rays - The rays, traced to the plate's surface, that you want to analyze
        
        Grating objects need the following Parameters:
        thickness - The thickness of the wires
        sep - The distance between the centers of adjacent wires
        
        Outputs:
        A trutharray, a photon's value is True if it collides with a wire
        
        Notes:
        - Since the wires run in the self.surface() direction, only the 
            x-position of the photons is considered
        - The "middle" wire is centered at x=0. Thus, a photon with x=0 will 
            always be removed.
        
        Example:
        >> c = CollimatorPlate(x=0,y=0,z=0,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0)
        >> c.collfunc = CollimatorPlate.wires
        >> c.thickness = .1
        >> c.sep = 1
        >> c.trace()
        '''
        if (type(self.thickness) != u.quantity.Quantity or 
        type(self.sep) != u.quantity.Quantity):
            raise ValueError('self.thickness and self.sep must be astropy units of length')
        self.sep = self.sep.to(u.mm)
        self.thickness = self.thickness.to(u.mm)
        
        x,y = self.getPosns(rays)
        
        # Condense the x-positions so they all exist on 0 < x < sep
        x = np.abs(x) % self.sep.value
        
        # Return True if the photon hits the wire to its left or right
        return np.logical_or((x < self.thickness.value/2),(x > self.sep.value - (self.thickness.value/2)))
        
    
    ## Removing Rays
    # Removing Rays is what a collimator does best
    
    def hit(self, rays):
        '''
        Function hit:
        This function simply return those rays that would be removed by this 
            collimator plate without actually removing them
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        
        Outputs:
        tarray - A trutharray, containing True if the photon has hit the plate 
            (and should be removed), contains False if the photon passes through 
            the plate (and should not be removed)
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up 
            to the user.
        - The function can use both a rectangle (from self.l and self.w) and a 
            collision function
        '''
        collfuncarray = np.zeros(len(rays))
        rectarray = np.zeros(len(rays))
        # Call the collision function (if it exists)
        if self.collisionfunction is not None:
            collfuncarray = self.collisionfunction(self,rays)
        # Check if rays miss the plate (if length and width are defined)
        if (self.l is not None) and (self.w is not None):
            x,y = self.getPosns(rays)
            rectarray = np.logical_or(np.abs(x) > self.w.value/2, np.abs(y) > self.l.value/2)
        # Return an array showing if photons should be removed by 
        # collisionfunction or missed the plate altogether
        return np.logical_or(collfuncarray,rectarray)
    
    def removemissed(self,rays,considerweights):
        '''
        Function removemissed:
        Removes the rays which have missed the collimator.
        This function is slightly different than Grating.removemissed.
        The Grating version will handle either the Grating's collision function 
            or the length and width, but not both.
        CollimatorPlate.removemissed can handle both a collision function and a 
            length/width.
        This is done to make the collision function simpler. Collimators will 
        already have complicated collision functions, so having removemissed 
        automatically handle length and width makes it much easier.
        
        Inputs:
        rays - a Rays Object which has been traced to this CollimatorPlate
        considerweights - Should be true if the photons are weighted
        
        
        Output:
        Up to two tuples containing how many rays were passed in and how many 
            survived. One tuple will show how many actually fell onto the 
            rectangle of the plate (if length and width are not None) The other 
            tuple will show how many photons survived the collisionfunction (if 
            it is not None)
        
        Notes:
        -self.collisionfunction should return True if the photon collides with 
            the collimator plate and needs to be removed
        '''
        # Store the length of the incoming rays
        l = rays.length(considerweights)
        
        # Initialize the efficiencies that we will eventually return
        # Any time we need to remove photons we will append to this list, so the
        # user can see specifically how rays were removed
        effs = []
        
        # get the positions of the rays:
        x,y = self.getPosns(rays)
        
        # Check if the photons missed the plate (if length and width are defined)
        if (self.l is not None) and (self.w is not None):
            tarray = np.logical_or(np.abs(x) > self.w.value/2, np.abs(y) > self.l.value/2)
            rays.remove(tarray)
            
            # Store the tuple we will eventually return
            effs.append(("Missed Collimator",l,rays.length(considerweights)))
        
        # Get the new length of the rays
        l = rays.length(considerweights)
        
        # apply the collision function if it exists
        if (self.collisionfunction is not None):
            tarray = self.collisionfunction(self,rays)
            
            rays.remove(tarray)
            
            # Store the tuple we will eventually return
            effs.append(("Eliminated by Collimator",l,rays.length(considerweights)))
        
        return effs
    
    
    ## Trace Function:
    
    def trace(self,rays,considerweights=False):
        '''
        Function trace:
        Traces rays to this Collimator Plate and removes photons as necessary. 
        This is a function that requires no input from the user and thus will be 
        called by the Instrument Object.
        
        Inputs:
        rays - The rays you want to trace to this plate
        
        Outputs:
        The efficiency, which tracks how many photons were removed by this
        plate
        '''
        self.trace_to_surf(rays)
        return self.removemissed(rays,considerweights)