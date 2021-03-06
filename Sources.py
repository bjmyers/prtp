import numpy as np
from prtp.Rays import Rays
import astropy.units as u

class Source:
    '''
    class Source:
    This class is what is initially passed to an Instrument object. This class
        contains information about how to generate a Rays object for use in
        an instrument
    '''
    
    def __init__(self, num = 0, wave = None, order = None):
        '''
        Function __init__:
        
        Inputs:
        num - The number of photons you want to generate
        wave - The wavelength of the photons. Can be a float or an array that
            has length num. Must be an astropy Quantity
        order - The order of the photons. Can be an integer or an array that
            has length num
        '''
        if (wave is None):
            self.wave = wave
        else:
            self.wave = wave.to(u.nm,equivalencies=u.spectral())
        self.num = num
        self.order = order
        self.rays = None
    
    @classmethod
    def sourceFromRays(cls, rays):
        '''
        Function sourceFromRays:
        An alternate initialization method. Takes in a Rays object and saves it.
        
        Inputs:
        rays - The Rays object you would like to eventually output.
        
        Notes:
        - This initialization method should be used when transferring rays 
            from one Instrument to another.
        '''
        source = Source()
        source.rays = rays
        return source
    
    @u.quantity_input(dx=u.mm,dy=u.mm,dz=u.mm,dl=u.rad,dm=u.rad,dn=u.rad)
    def addMisalignment(self, dx = 0*u.mm, dy = 0*u.mm, dz = 0*u.mm, 
                              dl = 0*u.rad, dm = 0*u.rad, dn = 0*u.rad):
        '''
        Function addMisalignment:
        Initializes misalingments of the Rays
        
        Inputs:
        dx,dy,... - Give the amount that the rays should be transformed in each
            direction. They must all be astropy units
        
        Outputs:
        None
        '''
        self.dx = dx.to(u.mm)
        self.dy = dy.to(u.mm)
        self.dz = dz.to(u.mm)
        self.dl = dl.to(u.rad)
        self.dm = dm.to(u.rad)
        self.dn = dn.to(u.rad)
    
    def misalign(self,rays):
        '''
        Function Misalign:
        Executes the transformation to misalign the rays
        
        Inputs:
        rays - The rays to be misaligned
        
        Outputs:
        Nothing
        
        Notes:
        - This function uses a try-except block to see if the misalignment 
            values have been specified (if they have not bee initialized, trying
            to access them will yield an error, which means we don't have to
            misalign the rays at all.) But this method makes it difficult to
            determine if a different error is occurring here, consider changing
            how this function works
        '''
        try:
            # Add misalignments if they've been defined
            rays.transform(self.dx.value,self.dy.value,self.dz.value,
                           self.dl.value,self.dm.value,self.dn.value)
        except:
            # This executes if misalignments have not been defined
            pass
    
    def addParams(self,rays):
        '''
        Function addParams:
        Adds the Wavelength and Order Parameters to the Rays object
        
        Inputs:
        rays - The Rays object you want to add Params to
        
        Outputs:
        None
        
        Notes:
        - If no parameters have been defined, errors will be caused down the line
            and it may be difficult to recognize them. But I also want the user
            to be able to define a Source without waves and orders. Think about
            adding a warning to the user that no waves or orders have been 
            added
        '''
        
        # Simplified Variant
        if (self.wave is not None):
            if type(self.wave.value) == np.ndarray:
                rays.wave = self.wave.value
            else:
                rays.wave = np.ones(len(rays)) * self.wave.value
        if (self.order is not None):
            if type(self.order) == np.ndarray:
                rays.order = self.order
            else:
                rays.order = np.ones(len(rays)) * self.order
    
    def generateRays(self):
        '''
        Function generateRays:
        This function should only be called on Source objects (as opposed to
            Source subclass objects) when it has been initialized already
            with a Rays object.
        '''
        if self.rays is not None:
            return self.rays
        else:
            raise ValueError('This Source has not been initialized with a Rays object')

class Annulus(Source):
    
    @u.quantity_input(rin=u.mm,rout=u.mm,z=u.mm)
    def __init__(self,num=1000,rin=1*u.mm,rout=2*u.mm,zhat=-1.,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        
        self.rin = rin.to(u.mm)
        self.rout = rout.to(u.mm)
        self.zhat = zhat
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.annulus(self.rin.value,self.rout.value,self.num,self.zhat)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays


class CircularBeam(Source):
    
    @u.quantity_input(rad=u.mm,z=u.mm)
    def __init__(self,num=1000,rad=1*u.mm,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.rad = rad.to(u.mm)
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.circularbeam(self.rad.value,self.num)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays


class ConvergingBeam(Source):
    
    @u.quantity_input(zset=u.mm,rin=u.mm,rout=u.mm,tmin=u.rad,tmax=u.rad,lscat=u.arcsec,z=u.mm)
    def __init__(self,num=1000,zset=0*u.mm,rin=0*u.mm,rout=1*u.mm,
    tmin=0*u.rad,tmax=1*u.rad,lscat=0*u.arcsec,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        
        self.zset = zset.to(u.mm)
        self.rin = rin.to(u.mm)
        self.rout = rout.to(u.mm)
        self.tmin = tmin.to(u.rad)
        self.tmax = tmax.to(u.rad)
        self.lscat = lscat.to(u.arcsec)
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.convergingbeam(self.zset.value,self.rin.value,
        self.rout.value,self.tmin.value,self.tmax.value,self.num,
        self.lscat.value)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays


class ConvergingBeam2(Source):
    
    @u.quantity_input(zset=u.mm,xmin=u.mm,xmax=u.mm,ymin=u.mm,ymax=u.mm,lscat=u.arcsec,z=u.mm)
    def __init__(self,num=1000,zset=0*u.mm,xmin=0*u.mm,xmax=1*u.mm,
    ymin=0*u.mm,ymax=1*u.mm,lscat=0*u.arcsec,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        
        self.zset = zset.to(u.mm)
        self.xmin = xmin.to(u.mm)
        self.xmax = xmax.to(u.mm)
        self.ymin = ymin.to(u.mm)
        self.ymax = ymax.to(u.mm)
        self.lscat = lscat.to(u.arcsec)
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.convergingbeam2(self.zset.value,self.xmin.value,
        self.xmax.value,self.ymin.value,self.ymax.value,self.num,
        self.lscat.value)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays


class PointSource(Source):
    
    @u.quantity_input(ang=u.rad,z=u.mm)
    def __init__(self,num=1000,ang=.1*u.rad,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.ang = ang.to(u.rad)
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.pointsource(self.ang.value,self.num)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays


class RectArray(Source):
    
    @u.quantity_input(xsize=u.mm,ysize=u.mm,z=u.mm)
    def __init__(self,num=100,xsize=1*u.mm,ysize=1*u.mm,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.xsize = xsize.to(u.mm)
        self.ysize = ysize.to(u.mm)
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.rectarray(self.xsize.value,self.ysize.value,self.num)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays


class RectBeam(Source):
    
    @u.quantity_input(xhalfwidth=u.mm,yhalfwidth=u.mm,z=u.mm)
    def __init__(self,num=1000,xhalfwidth=1*u.mm,yhalfwidth=1*u.mm,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.xhalfwidth = xhalfwidth.to(u.mm)
        self.yhalfwidth = yhalfwidth.to(u.mm)
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.rectbeam(self.xhalfwidth.value,self.yhalfwidth.value,self.num)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays


class Subannulus(Source):
    
    @u.quantity_input(rin=u.mm,rout=u.mm,dphi=u.rad,z=u.mm)
    def __init__(self,num=1000,rin=0*u.mm,rout=1*u.mm,dphi=1*u.rad,zhat=1.,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.rin = rin.to(u.mm)
        self.rout = rout.to(u.mm)
        self.dphi = dphi.to(u.rad)
        self.zhat = zhat
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.subannulus(self.rin.value,self.rout.value,self.dphi.value,
        self.num,self.zhat)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays


class Xslit(Source):
    
    @u.quantity_input(xin=u.mm,xout=u.mm,z=u.mm)
    def __init__(self,num=1000,xin=0*u.mm,xout=1*u.mm,zhat=-1.,z=0*u.mm,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.xin = xin.to(u.mm)
        self.xout = xout.to(u.mm)
        self.zhat = zhat
        self.z = z.to(u.mm)
    
    def generateRays(self):
        rays = Rays.xslit(self.xin.value,self.xout.value,self.num,self.zhat)
        rays.translate(dz=self.z.value)
        self.addParams(rays)
        self.misalign(rays)
        return rays