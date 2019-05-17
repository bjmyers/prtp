import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
from prtp.FlatComponent import FlatComponent
import prtp.transformationsf as trans

class Grating(FlatComponent):

    
    ## Initialization Functions
    
    def __init__(self,x=0,y=0,z=0,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0, l=1,w=1,pfunc=None,collfunc=None):
        '''
        Initializes a Grating Object, requires the following arguments:
        
        x,y,z - The position of a point along the grating
        nx,ny,nz - The components of a vector normal to the grating surface
        sx,sy,sz - The surface vector of this grating, the mean groove direction
        l - The length of the grating, that is, its extent in the s direction
        w - The width of the grating, that is, its extent in the sxn direction
        pfunc - A user-defined function to determine the groove period experienced by each photon
        collfunc - A user-defined function to determine if Rays have missed the Component
        '''
        FlatComponent.__init__(self,x,y,z,nx,ny,nz,fx,fy,fz, collfunc=collfunc)
        self.l = l
        self.w = w
        self.periodfunction = pfunc
    
    
    def copy(self):
        '''
        Function copy:
        Returns a copy of this Grating Object
        
        Inputs:
        None
        
        Outputs:
        An identical Grating Object with the same attributes
        '''
        return Grating(self.x,self.y,self.z,
        self.nx,self.ny,self.nz,
        self.fx,self.fy,self.fz,
        self.l,self.w,self.periodfunction,self.collisionfunction)

    
    ## Analyzing Traced Rays
    
    def hit(self,rays):
        '''
        Function hit:
        Given Rays that have been traced to the Grating, returns a trutharray detailing which photons have fallen within the rectangular region of this grating.
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        
        Outputs:
        tarray - A trutharray, containing True if the photon has fallen within the grating's dimensions or False if the photon missed the grating
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        - The function use a rectangle defined by length and width (self.l and self.w) or a user-defined collision function (self.collfunc). The user-defined function will be performed preferentially.
        '''
        if self.periodfunction is not None:
            return self.periodfunction(rays)
        elif (self.l is not None) and (self.w is not None):
            x,y = self.getPosns(rays)
            return np.logical_and(np.abs(x) < self.w/2, np.abs(y) < self.l/2)
        else:
            raise NotImplementedError('This Grating needs a collision function to be implemented or a length and width if you wish to use a simple rectangular model')
    
    def removemissed(self,rays):
        '''
        Function removemissed:
        Given Rays that have been traced to the Grating, returns the photons which have hit the Grating.
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        
        Outputs:
        A tuple containing the original number of rays and the final number of rays
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        - Make sure the length and the width have been defined when the Grating was defined
        '''
        l = len(rays)
        # Find rays which have not hit
        tarray = np.logical_not(self.hit(rays))
        rays.remove(tarray)
        return (l,l-len(rays))
    
    
    ## Sample Period Functions
    # Some already defined functions in case you don't want to define your own
    def linear(self, rays, period, y=0, focus=100):
        '''
        Function linear:
        A sample period function that simulates linearly converging grooves.
        
        Inputs:
        rays - The rays that have been traced to the grating.
        period - The period known at one position
        y - The position at which the period is period. By default it is the center of the grating
        focus - The position of the focus, at this point the period is 0
        
        Outputs:
        An array of periods for each input photon
        
        Notes:
        The period assigned by this function only depends on a photon's y-position.
        '''
        raysx,raysy = self.getPosns(rays)
        return period - (period)/(focus-y) * (raysy - y)

        
    
    ## Grating Functions:
    # Handle reflecting off of the gratings
    
    def getPeriod(self, rays, **kwargs):
        '''
        Function getPeriod:
        Calls the user-defined periodfunction. This is just an easy place to store and call your function to obtain the period for each photon
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        **kwargs - Handles any other arguments a user-defined function might require. For example, a function might require a grating period at a single point, this argument would be named in the periodfunction and passed as an argument here. 
        
        Outputs:
        An array of periods, with each grating period corresponding to a photon
        
        Notes:
        - User-Defined period functions must take in Rays objects and return arrays of grating periods experienced by each photon. Be sure to 
        - The output can be passed into self.grat or self.radgrat as the dor dpermm argument, respectively
        - If no period function has been defined, this function will return None
        '''
        if self.periodfunction is not None:
            return self.periodfunction(rays,**kwargs)
    
    def grat(self, rays, d=160, order=0, wave=100, autoreflect=False):
        '''
        Function grat:
        Given Rays that have been traced to the Grating, reflects the Rays off of a parallel grating
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        d - The period of the grooves, can be a float or array
        order - The order of the photon, can be a float or array
        wave - the wavelength of the photon, can be a float or array
        autoreflect - Boolean. If true, the photons will automatically be reflected off of the surface
        
        Outputs:
        rays - The updated Rays object that has been reflected off of the Grating
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        - Photons must be reflected after this function is called or they will not behave properly
        '''
        vel = [rays.l,rays.m,rays.n]
        # Define normal vector, focus vector, and their cross product, stacked on themselves so they can be dotted into the rays
        nor = self.Normal()
        s = self.Surface()
        sxn = np.cross(s,n)
        norm = np.tile(n,(len(rays),1)).transpose()
        surf = np.tile(s,(len(rays),1)).transpose()
        crossproduct = np.tile(np.cross(s,n),(len(rays),1)).transpose()
        
        # Define velocity components aligned with Grating (if it were in the xy-plane) using dot products
        gratl = (vel * crossproduct).sum(0)
        gratm = (vel * surf).sum(0)
        gratn = (vel * norm).sum(0)
        
        # This gives us the new velocity in grating components, we need to convert it to xyz components
        gratl,gratn,gratm = trans.grat(gratl,gratm,gratn,d,order,wave)
        
        # Assumes normal and surface vectors are normalized, they should be unless the user has been messing stuff up
        l = gratl*sxn[0] + gratm*s[0] + gratn*nor[0]
        m = gratl*sxn[1] + gratm*s[1] + gratn*nor[1]
        n = gratl*sxn[2] + gratm*s[2] + gratn*nor[2]
        
        rays.set(l=l,m=m,n=n)
        
        if (autoreflect):
            self.reflect(rays)
    
    def radgrat(self,rays,dpermm=160,order=0,wave=100,autoreflect=False):
        '''
        Function radgrat:
        Given Rays that have been traced to the Grating, reflects the Rays off of a radial grating
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        dpermm - The period of the grooves, can be a float or array
        order - The order of the photon, can be a float or array
        wave - the wavelength of the photon, can be a float or array
        autoreflect - Boolean. If true, the photons will automatically be reflected off of the surface
        
        Outputs:
        rays - The updated Rays object that has been reflected off of the Grating
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        - Photons must be reflected after this function is called or they will not behave properly
        '''
        vel = [rays.l,rays.m,rays.n]
        # Define normal vector, focus vector, and their cross product, stacked on themselves so they can be dotted into the rays
        nor = self.Normal()
        s = self.Surface()
        sxn = np.cross(s,nor)
        norm = np.tile(nor,(len(rays),1)).transpose()
        surf = np.tile(s,(len(rays),1)).transpose()
        crossproduct = np.tile(np.cross(s,nor),(len(rays),1)).transpose()
        
        # Define velocity components aligned with Grating (if it were in the xy-plane) using dot products
        gratl = (vel * crossproduct).sum(0) # Vel aligned with sxn vector
        gratm = (vel * surf).sum(0)         # Vel aligned with surf vector
        gratn = (vel * norm).sum(0)         # Vel aligned with norm vector
        
        # This gives us the new velocity in grating components, we need to convert it to xyz components
        x,y,l,m,n = trans.radgrat(rays.x,rays.y,gratl,gratm,gratn,dpermm,order,wave)
        
        # Assumes normal and surface vectors are normalized, they should be unless the user has been messing stuff up
        l = gratl*sxn[0] + gratm*s[0] + gratn*nor[0]
        m = gratl*sxn[1] + gratm*s[1] + gratn*nor[1]
        n = gratl*sxn[2] + gratm*s[2] + gratn*nor[2]
        
        rays.set(l=l,m=m,n=n)
        
        if (autoreflect):
            self.reflect(rays)
    
    def reflect(self,rays):
        '''
        Function reflect:
        Given Rays that have been traced to the Grating, reflectes the Rays off of the surface. This function is different than grat or radgrat, and needs to be called if you want the photons to be accurately reflected off the Grating
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        dpermm - The period of the grooves, can be a float or array
        order - The order of the photon, can be a float or array
        wave - the wavelength of the photon, can be a float or array
        
        Outputs:
        rays - The updated Rays object that has been reflected off of the Grating
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        - Photons must be reflected after this function is called or they will not behave properly
        '''
        length = len(rays)
        rays.ux = np.ones(length) * self.nx
        rays.uy = np.ones(length) * self.ny
        rays.uz = np.ones(length) * self.nz
        rays.reflect()














