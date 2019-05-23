import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
from prtp.FlatComponent import FlatComponent
import prtp.transformationsf as trans

class Grating(FlatComponent):

    
    ## Initialization Functions
    
    def __init__(self,x=0,y=0,z=0,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0,
    l=None,w=None,pfunc=None,collfunc=None,radial=True,d=160,fdist=None):
        '''
        Initializes a Grating Object, requires the following arguments:
        
        x,y,z - The position of a point along the grating
        nx,ny,nz - The components of a vector normal to the grating surface
        sx,sy,sz - The surface vector of this grating, the mean groove direction
        l - The length of the grating, that is, its extent in the s direction
        w - The width of the grating, that is, its extent in the sxn direction
        pfunc - A user-defined function to determine the groove period 
            experienced by each photon. If the pfunc takes arguments other than
            just rays, these need to be defined as Grating parameters. See the
            documentation for Grating.linear for an example
        collfunc - A user-defined function to determine if Rays have missed the 
            Grating
        radial - If True, this grating will behave as a radial grating.
            If False, this grating will behave as a paralle grating
        d - The period of the Grating, a pfunc can be used for more sophisticated
            periods. For a radial grating, it is the period at the center of the
            Grating
        fdist - The distance to the focus of the grooves. It should be measured
            from the center of the grating to the focus along the s-direction
            This value is important for radial Gratings where we need the 
            distance of every photon from the groove focus.
        '''
        FlatComponent.__init__(self,x,y,z,nx,ny,nz,sx,sy,sz, collfunc=collfunc)
        self.l = l
        self.w = w
        self.periodfunction = pfunc
        self.collfunc = collfunc
        self.radial = radial
        self.d = d
        self.fdist = fdist
    
    
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
        self.sx,self.sy,self.sz,
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
        if self.collfunc is not None:
            return self.collfunc(rays)
        elif (self.l is not None) and (self.w is not None):
            x,y = self.getPosns(rays)
            return np.logical_and(np.abs(x) < self.w/2, np.abs(y) < self.l/2)
        else:
            raise NotImplementedError('This Grating needs a collision function to be implemented or a length and width if you wish to use a simple rectangular model')
    
    def removemissed(self,rays,considerweights=False):
        '''
        Function removemissed:
        Given Rays that have been traced to the Grating, returns the photons which have hit the Grating.
        
        Inputs:
        rays - A Rays object that has been traced to the Grating Plane
        considerweights - Should be used if the photons are weighted
        
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
        return ("Missed Grating",l,rays.length(considerweights=considerweights))
    
    
    ## Sample Period Functions
    # Some already defined functions in case you don't want to define your own
    def linear(self, rays):
        '''
        Function linear:
        A sample period function that simulates linearly converging grooves.
        
        Inputs:
        rays - The rays that have been traced to the grating.
        
        Grating objects need the following Parameters:
        d - The period known at one position
        y0 - The position in the y-direction where the period is d
        focus - The position of the focus, at this point the period is 0
        
        Outputs:
        An array of periods for each input photon
        
        Notes:
        The period assigned by this function only depends on a photon's y-position.
        
        Example of Use:
        >> g = Grating(x=0,y=0,z=0,nx=0,ny=0,nz=1,sx=0,sy=1,sz=0,d=160)
        >> g.y0 = 0
        >> g.focus = 1e8
        >> g.trace_to_surf(rays)
        >> g.radgrat(rays,d=Grating.linear(rays))
        '''
        raysx,raysy = self.getPosns(rays)
        return self.d - (self.d)/(self.focus-self.y0) * (raysy - self.y0)
        
    
    ## Grating Functions:
    # Handle reflecting off of the gratings
    
    def grat(self, rays, order=0, wave=1., autoreflect=False):
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
        
        # Find the groove period each photon sees
        if self.periodfunction is None:
            d = self.d
        else:
            d = self.periodfunction(rays)
        
        # This gives us the new velocity in grating components, we need to convert it to xyz components
        gratl,gratn,gratm = trans.grat(gratl,gratm,gratn,d,order,wave)
        
        # Assumes normal and surface vectors are normalized, they should be unless the user has been messing stuff up
        l = gratl*sxn[0] + gratm*s[0] + gratn*nor[0]
        m = gratl*sxn[1] + gratm*s[1] + gratn*nor[1]
        n = gratl*sxn[2] + gratm*s[2] + gratn*nor[2]
        
        rays.set(l=l,m=m,n=n)
        
        if (autoreflect):
            self.reflect(rays)
    
    def radgrat(self,rays,order=0,wave=1.,autoreflect=False):
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
        
        # Find the groove period each photon sees
        if self.periodfunction is None:
            d = self.d
        else:
            d = self.periodfunction(rays)
        
        
        x,y = self.getPosns(rays)
        y -= self.fdist
        d /= self.fdist
        
        # This gives us the new velocity in grating components, we need to convert it to xyz components
        x,y,l,m,n = trans.radgrat(x,y,gratl,gratm,gratn,d,order,wave)
        
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
    
    
    def trace(self,rays,considerweights=False):
        '''
        Function trace:
        Traces rays to and reflects them off of this Grating. This is a function
        that requires no input from the user and thus will be called by the
        Instrument Object
        
        Inputs:
        rays - The rays you want to trace to this grating
        
        Outputs:
        eff - The efficiency, currently only tracks how many photons missed the
        grating. In the future will consider absorption from the material too.
        '''
        self.trace_to_surf(rays)
        eff1 = self.removemissed(rays,considerweights=considerweights)
        if self.radgrat:
            self.radgrat(rays,order=rays.getParam('Order'),wave=rays.getParam('Wavelength'),autoreflect=True)
            # Need to check l for nans b/c rays.x is not modified
            rays.remove(np.isnan(rays.l))
        else:
            self.grat(rays,order=rays.getParam('Order'),wave=rays.getParam('Wavelength'),autoreflect=True)
            rays.remove(np.isnan(rays.l))
            
        eff2 = ('Failed to Reflect off Grating', eff1[2],rays.length(considerweights))
        return [eff1,eff2]














