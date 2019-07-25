import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
from prtp.FlatComponent import FlatComponent
import prtp.transformationsf as trans
import astropy.units as u
from scipy.integrate import quad

class Detector(FlatComponent):
    
    ## Initialization Functions:
    
    @u.quantity_input(x=u.mm,y=u.mm,z=u.mm,l=u.mm,w=u.mm,fieldfree=u.um,depletion=u.um)
    def __init__(self,x=0*u.mm,y=0*u.mm,z=0*u.mm,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0,
        q=1.,l=1*u.mm,w=1*u.mm,xpix=10,ypix=10,fieldfree=15*u.um,depletion=3*u.um,
        considersplits=False, steps=50):
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
        fieldfree - The depth of the fieldfree region of the detector, default
            units are in microns, must be in units of length.
        depletion - The depth of the depletion region of the detector, is not
            currently used. Must be in units of length, default units are 
            microns.
        considersplits - A boolean to determine if split events should be 
            considered. Should be True if you want to simulate split events,
            note that this will slow down the view() function significantly.
        steps - If considersplits is True, the split events will be calculated
            using this number of steps. Lower number of steps will be less 
            accurate but will run faster. For small splits (taking up only a
            few pixels), few steps is acceptable.
        
        Notes:
        -Unlike other FlatComponents, Detector needs a length and width, those 
            two parameters cannot be None
        '''
        FlatComponent.__init__(self,x,y,z,nx,ny,nz,sx,sy,sz)
        self.q = q
        self.l = l.to(u.mm)
        self.w = w.to(u.mm)
        self.xpix = xpix
        self.ypix = ypix
        self.fieldfree = fieldfree
        self.depletion = depletion
        self.considersplits = considersplits
        self.steps = steps
        
        # Initialize the pixels array
        self.pixels = np.zeros((xpix,ypix))
        
        # This param tells what type of noise is added in view()
        self.noise = None
        
        # This param stores how many custom noise frames have been added
        self.numframes= 0
    
    
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
        tarray - A trutharray, containing True if the photon has hit the 
            detector, contains False if the photon missed the detector
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        '''
        x,y = self.getPosns(rays)
        return np.logical_and(np.abs(x) < self.w.value/2, np.abs(y) < self.l.value/2)
        
    
    def removemissed(self,rays,considerweights=False,eliminate='remove'):
        '''
        Function removemissed:
        Removes the rays which have missed the Detector. Also removes some 
            photons according to the quantum efficiency (e.g: if self.q = .1, 
            10% of the photons that hit will be removed)
        
        Inputs:
        rays - a Rays Object which has been traced to this CollimatorPlate
        considerweights - should True if the photons are weighted
        eliminate - If 'remove', photons that miss will be removed. Otherwise,
            missed photons will be replaced with NaNs in the x-position
        
        Output:
        Two tuples describing the efficiency of the detector:
            The first one contains (num_input_rays,num_hit_rays)
            The second one contains (num_hit_rays,num_detected_rays)
        '''
        
        l = rays.length(considerweights)
        # Find rays which have not hit
        tarray = np.logical_not(self.hit(rays))
        
        if eliminate == 'remove':
            rays.remove(tarray)
        else:
            rays.x[tarray] = np.nan
        
        t1 = ("Missed Detector",l,rays.length(considerweights))
        l = rays.length(considerweights)
        
        rays.probRemove(self.q)
        t2 = ("Eliminated by QE of Detector",l,rays.length(considerweights))
        
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
        Tells the detector that Gaussian Noise must be added when the photons 
            are viewed
        
        Inputs:
        mean - The mean of the Gaussian Distribution
        std - The Standard Deviation of the Gaussian Distribution
        
        Outputs:
        None
        '''
        self.noise = 'gaussian'
        self.mean = mean
        self.std = std
    
    def addNoise(self,noise):
        '''
        Adds a custom noise profile to this detector.
        
        Inputs:
        noise - A sample frame containing only noise, could be a dark frame,
            bias frame, flat field, etc: It must be a 2D array with the same
            shape as the pixel array. The values of the noise input are saved,
            and if a noise frame has been added previously, the values are 
            averaged.
        
        Outputs:
        None
        
        Notes:
        - When view() is called and self.noise is 'custom', noise will be added
            to the pixel array in a poisson distribution with lambda equal to 
            the values given by noisevalues.
        - Noise Frames must be added consecutively, called addGaussianNoise()
            will delete any noise frames that have been previously added. For
            example, if you want to average 10 noise frames, you need ten
            calls to addNoise without calling addGaussianNoise at all.
        '''
        # Refresh numframes if noise was not custom before
        if self.noise != 'custom':
            self.numframes = 0
            self.noisevalues = np.zeros((self.xpix,self.ypix))
        
        self.noise = 'custom'
        
        # Using previous means, calculate the new mean when the inputted noise
        # is added to the dataset
        self.noisevalues *= self.numframes
        self.noisevalues += noise
        self.noisevalues /= (self.numframes + 1)
        
        self.numframes += 1
    
    
    ## Viewing Functions:
    # See what the output of the detector will look like
        
    
    def view(self, rays=None):
        '''
        Function view:
        Shows what the readout of the Detector should look like
        
        Inputs:
        rays - A Rays object that has been traced to the Detector and has had 
            the missed rays removed
            if Rays is None, this function will simply return the current 
                pixel array
        
        Output:
        temppixels - A pixel array of the Detector readout, can be viewed 
            easily with plt.imshow()
        
        Notes:
        - The original self.pixels is unmodified by this function
        - This function will consider split events if self.considersplits is
            True
        '''
        if rays is None:
            return self.pixels
        
        if self.noise == 'gaussian':
            self.pixels += np.random.normal(loc=self.mean,scale=self.std,
                                            size=self.pixels.shape)
        if self.noise == 'custom':
            self.pixels += np.random.poisson(lam=self.noisevalues,
                                            size=self.pixels.shape)
        
        if (self.considersplits):
            self.pixels += self.split(rays,self.steps)
        
        else:
        
            x,y = self.getPosns(rays)
            
            # Move axes so that the origin is in the bottom-left corner
            x += self.l.value/2
            y += self.w.value/2
            
            # Find out the length and width of each pixel
            xperpix = self.l.value / self.xpix
            yperpix = self.w.value / self.ypix
            
            xposns = (x // xperpix).astype(int)
            yposns = (y // yperpix).astype(int)
            
            for i in range(len(xposns)):
                self.pixels[xposns[i],yposns[i]] += (rays.wave[i]*u.nm).to(u.eV,equivalencies=u.spectral()).value/3.65
        
        return self.pixels
    
    
    
    def split(self, rays, steps = 50):
        '''
        Function split:
        Calculates how the charge cloud produced by each photon spill out into
            adjacent pixels. The algorithm works by considering small slices
            of the charge cloud and figuring out when the slice moves into 
            adjacent pixels. Using the area of the slice and a Gaussian integral,
            each slice will deposit charge into nearby pixels.
        
        Inputs:
        rays - The photons you want to simulate
        steps - The number of slices you want to simulate. A smaller number is
            faster but less accurate. Accuracy is most important for large 
            charge clouds.
        
        Outputs:
        arr - An array the same size as the Detector's pixel array, containing
            the charge cloud binned into pixels. This array is added to the 
            noise array by the view() function.
        '''
        # Get the x and y positions of the photons
        xs = rays.x
        ys = rays.y
        
        # Offset the positions so the origin is in the corner
        # This is also done for the noise matrix
        xs += self.l.value/2
        ys += self.w.value/2
        
        # Get the energies of the photons in eV
        energies = (rays.wave*u.nm).to(u.eV,equivalencies=u.spectral()).value
        
        # Find the number of electrons created by dividing the energy by the
        # quantum yield of Silicon
        counts = energies / 3.65
        
        # Find the length and width of each pixel
        pixlen = self.l.value / self.xpix
        pixwid = self.w.value / self.ypix
        
        # Fine the attenuation length of silicon at this energy in um
        # Found using a 2nd order polynomial approximation to the attenuation
        # length as a function of energy
        attenlens = 5.86902987e-06*energies**2 + -4.128074e-03*energies + 9.79687932e-01
        
        # Find how far the photons penetrate into the Silicon using the 
        # attenuation length and the field free depth
        d = (self.fieldfree.to(u.um) - (attenlens*u.um)).to(u.mm).value
        
        # Charge cloud is a Gaussian with FWHM=2d, find the standard deviation
        # of the distribution
        sigma = d / (np.sqrt(2*np.log(2)))
        
        # Initialize the list of slices we will be using in the future
        thetas = np.linspace(0.0001,2*np.pi,steps)
        
        # Find the angle that once slice subtends
        arc = np.mean(np.diff(thetas) / (2*np.pi))
        
        # Initialize a blank initial array
        arr = np.zeros((self.xpix,self.ypix))
    
        for i in range(len(thetas)):
            
            theta = thetas[i]
            costheta = np.cos(theta)
            sintheta = np.sin(theta)
            
            # Find the coordinates of the pixel that each photon starts at
            xpos = (xs/pixlen).astype(int)
            ypos = (ys/pixwid).astype(int)
            
            # Find the offset of each photon's initial position from the center
            # of the pixel
            dx = pixlen/2 - (xs % pixlen)
            dy = pixwid/2 - (ys % pixwid)
            
            # We need to keep track of how many vertical and horizontal steps
            # that each photon has undergone so far
            num_horiz_steps = np.zeros(len(xs))
            num_vert_steps = np.zeros(len(xs))
            
            # This will keep track of the distance between the cloud's center 
            # and the previous pixel boundary
            prev_r = np.zeros(len(xs))
            
            # X and Y step are either 1 or -1 that depends on the direction of theta
            ystep = int(costheta / np.abs(costheta))
            xstep = int(sintheta / np.abs(sintheta))
            
            # This will keep track of the distance between the cloud's center 
            # and the next pixel boundary
            r = np.zeros(len(xs))
            
            # Keeps track of the photons which are done
            done = np.zeros(len(xs)).astype(bool)
            
            while True:
                
                # calculate the two potential step sizes
                next_horiz_step = (((2*num_horiz_steps + 1) * pixwid * ystep / 2) - dx) / np.abs(costheta)
                next_vert_step = (((2*num_vert_steps + 1) * pixlen * xstep / 2) - dy) / np.abs(sintheta)
                
                # Find which photons are done
                newdone = np.logical_or(np.logical_or(xpos < 0, xpos >= self.xpix),
                                    np.logical_or(ypos < 0, ypos >= self.ypix))
                done = np.logical_or(newdone,done)
                
                # Find the photons undergoing their last step
                # ------------------------------------------
                laststep = np.logical_and(np.abs(next_horiz_step) > d, np.abs(next_vert_step) > d)
                
                # Block the photons which are done
                laststep = np.logical_and(laststep,np.logical_not(done))
                
                # Update the positions of the laststep photons
                r[laststep] = d[laststep]
                
                # This is the indices of the photons undergoing their last step
                ind = np.where(laststep)[0]
                
                # Iterate through the last indices to find how many counts need to 
                # be added to this pixel
                for i in ind:
                    arr[xpos[i],ypos[i]] += quad(lambda x: weightedgauss(x), prev_r[i]/sigma[i], r[i]/sigma[i])[0] * energies[i] / steps
                
                # Remember that these photons do not need to be iterated anymore
                done[laststep] = 1
                
                # Stop if all photons are done
                if done.all():
                    break
                
                # Find the photons which need to take a horizontal step
                # --------------------------------------------------------
                needhoriz = np.abs(next_horiz_step) <= np.abs(next_vert_step)
                
                # Block the photons which are done
                needhoriz = np.logical_and(needhoriz,np.logical_not(done))
                
                # update the positions of the photons
                r[needhoriz] = next_horiz_step[needhoriz]
                
                # This is the indices of the photons undergoing horiz steps
                ind = np.where(needhoriz)[0]
                
                # Iterate through the horiz indices to find how many counts need to 
                # be added to this pixel
                for i in ind:
                    arr[xpos[i],ypos[i]] += quad(lambda x: weightedgauss(x), prev_r[i]/sigma[i], r[i]/sigma[i])[0] * energies[i] / steps
                
                # Update ypositions
                ypos[needhoriz] += ystep
                num_horiz_steps[needhoriz] += 1
                
                # Find the photons which need to take a vertical step
                # --------------------------------------------------------
                needvert = np.abs(next_vert_step) < np.abs(next_horiz_step)
                
                # Block the photons which are done
                needvert = np.logical_and(needvert,np.logical_not(done))
                
                # update the positions of the photons
                r[needvert] = next_vert_step[needvert]
                
                # This is the indices of the photons undergoing horiz steps
                ind = np.where(needvert)[0]
                
                # Iterate through the horiz indices to find how many counts need to 
                # be added to this pixel
                for i in ind:
                    arr[xpos[i],ypos[i]] += quad(lambda x: weightedgauss(x), prev_r[i]/sigma[i], r[i]/sigma[i])[0] * energies[i] / steps
                
                # Update ypositions
                xpos[needvert] += xstep
                num_vert_steps[needvert] += 1
                
                prev_r = r.copy()
        
        # Normalize Gaussian
        arr *= 5.013256549262
        
        return arr

    ## Trace Function:
    
    def trace(self,rays,considerweights=False,eliminate='remove'):
        '''
        Function trace:
        Traces rays to this Detector and removes photons as necessary. 
        This is a function that requires no input from the user and thus will be 
        called by the Instrument Object.
        
        Inputs:
        rays - The rays you want to trace to this detector
        considerweights - If true, any effect that probabilistically removes
            photons will instead affect their weights
        eliminate - If 'remove', photons that miss will be removed. Otherwise,
            missed photons will be replaced with NaNs in the x-position
        
        Outputs:
        The efficiency, which tracks how many photons were removed by this
        detector
        '''
        self.trace_to_surf(rays)
        return self.removemissed(rays,considerweights,eliminate)

@np.vectorize
def weightedgauss(r):
    '''
    Weighted Gauss is a function that needs to be integrated by the split() 
        function. If a Gaussian function is defined as G(r), the weighted
        gaussian is W(r) = r*G(r).
    '''
    # Assumes mean of 0 and a standard deviation of zero
    return r * np.exp(-.5 * r**2) / np.sqrt(2 * np.pi)


