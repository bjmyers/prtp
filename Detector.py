import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
from prtp.FlatComponent import FlatComponent
import prtp.transformationsf as trans
import astropy.units as u

class Detector(FlatComponent):
    
    ## Initialization Functions:
    
    def __init__(self,x=0,y=0,z=0,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0,q=1.,l=1,w=1,xpix=10,ypix=10):
        '''
        Initializes a Detector Object, requires the following arguments:
        
        x,y,z - The position of a point along the detector, must be astropy
            units of length
        nx,ny,nz - The components of a vector normal to the detector's surface
        sx,sy,sz - The surface vector of this plate, should be defined in a 
            direction that makes sense for the system (straight up, for example)
        q - The quantum efficiency of the detector, if q=.5, only 50% of the 
            rays that reach the detector will be observed
        l - The length of the detector, that is, its extent in the s direction,
            must be an astropy unit of length
        w - The width of the detector, that is, its extent in the sxn direction,
            must be an astropy unit of length
        xpix - The number of pixels along the length of the detector
        ypix - The number of pixels along the width of the detector
        
        Notes:
        -Unlike other FlatComponents, Detector needs a length and width, those 
            two parameters cannot be None
        '''
        FlatComponent.__init__(self,x,y,z,nx,ny,nz,sx,sy,sz)
        self.q = q
        if (type(l) != u.quantity.Quantity or type(w) != u.quantity.Quantity):
            raise ValueError('l and w must be astropy units of length')
        self.l = l.to(u.mm)
        self.w = w.to(u.mm)
        self.xpix = xpix
        self.ypix = ypix
        
        # Initialize the pixels array
        self.pixels = np.zeros((xpix,ypix))
    
    
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
        self.q,self.l,self.w,self.xpix,self.ypix)

    
    ## Removing Rays
    # Removing Rays is what a collimator does best
    
    def hit(self, rays):
        '''
        Function hit:
        This function returns the Rays that have hit the detector
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        
        Outputs:
        tarray - A trutharray, containing True if the photon has hit the detector, contains False if the photon missed the detector
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        '''
        x,y = self.getPosns(rays)
        return np.logical_and(np.abs(x) < self.w.value/2, np.abs(y) < self.l.value/2)
        
    
    def removemissed(self,rays,considerweights=False):
        '''
        Function removemissed:
        Removes the rays which have missed the Detector. Also removes some photons according to the quantum efficiency (e.g: if self.q = .1, 10% of the photons that hit will be removed)
        
        Inputs:
        rays - a Rays Object which has been traced to this CollimatorPlate
        considerweights - should True if the photons are weighted
        
        Output:
        Two tuples describing the efficiency of the detector:
            The first one contains (num_input_rays,num_hit_rays)
            The second one contains (num_hit_rays,num_detected_rays)
        '''
        
        l = rays.length(considerweights)
        # Find rays which have not hit
        tarray = np.logical_not(self.hit(rays))
        rays.remove(tarray)
        
        t1 = ("Missed Detector",l,rays.length(considerweights))
        l = rays.length(considerweights)
        
        rays.probRemove(self.q)
        t2 = ("Eliminated by QE of Detector",l,rays.length(considerweights=True))
        
        return [t1,t2]
    
    
    ## Reset:
    
    def reset(self):
        '''
        Function reset:
        Sets the value of all the pixels in the pixel array to 0
        
        Inputs:
        None
        
        Outputs:
        None
        '''
        self.pixels = np.zeros((self.xpix,self.ypix))
    
    
    ## Noise Functions:
    # add noise to the Detector
    
    def addGaussianNoise(self,mean=0,std=1):
        '''
        Function addGaussianNoise:
        Adds Gaussian noise to the pixel array
        
        Inputs:
        mean - The mean of the Gaussian Distribution
        std - The Standard Deviation of the Gaussian Distribution
        
        Outputs:
        None
        '''
        self.pixels += np.random.normal(loc=mean,scale=std,size=self.pixels.shape)
    
    
    ## Viewing Functions:
    # See what the output of the detector will look like
        
    
    def view(self, rays=None):
        '''
        Function view:
        Shows what the readout of the Detector should look like
        
        Inputs:
        rays - A Rays object that has been traced to the Detector and has had the missed rays removed
            if Rays is None, this function will simply return the current pixel array
        
        Output:
        temppixels - A pixel array of the Detector readout, can be viewed easily with plt.imshow()
        
        Notes:
        - The original self.pixels is unmodified by this function
        '''
        if rays is None:
            return self.pixels
        
        x,y = self.getPosns(rays)
        
        # Move axes so that the origin is in the bottom-left corner
        x += self.l.value/2
        y += self.w.value/2
        
        # Find out the length and width of each pixel
        xperpix = self.l.value / self.xpix
        yperpix = self.w.value / self.ypix
        
        xposns = (x // xperpix).astype(int)
        yposns = (y // yperpix).astype(int)
        
        temppixels = np.copy(self.pixels)
        
        # I'm sure there's a faster way to do this:
        for i in range(len(xposns)):
            temppixels[xposns[i]][yposns[i]] += 1
        
        return temppixels

    ## Trace Function:
    
    def trace(self,rays,considerweights=False):
        '''
        Function trace:
        Traces rays to this Detector and removes photons as necessary. 
        This is a function that requires no input from the user and thus will be 
        called by the Instrument Object.
        
        Inputs:
        rays - The rays you want to trace to this detector
        
        Outputs:
        The efficiency, which tracks how many photons were removed by this
        detector
        '''
        self.trace_to_surf(rays)
        return self.removemissed(rays,considerweights)
    
    
    
  



