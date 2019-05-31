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
            if (type(wave) != u.quantity.Quantity):
                raise ValueError('Wave must be an astropy Quantity')
            self.wave = wave.to(u.nm,equivalencies=u.spectral())
        self.num = num
        self.order = order
        self.rays = None
    
    def addMisalignment(self, dx = 0*u.mm, dy = 0*u.mm, dz = 0*u.mm, 
                              dl = 0*u.mm, dm = 0*u.mm, dn = 0*u.mm):
        '''
        Function addMisalignment:
        Initializes misalingments of the Rays
        
        Inputs:
        dx,dy,... - Give the amount that the rays should be transformed in each
            direction. They must all be astropy units
        
        Outputs:
        None
        '''
        if ((type(dx) != u.quantity.Quantity) or (type(dy) != u.quantity.Quantity)
        or (type(dz) != u.quantity.Quantity) or (type(dl) != u.quantity.Quantity)
        or (type(dm) != u.quantity.Quantity) or (type(dn) != u.quantity.Quantity)):
            raise ValueError('Any arguments of addMisalignment must be astropy units!')
        self.dx = dx.to(u.mm)
        self.dy = dy.to(u.mm)
        self.dz = dz.to(u.mm)
        self.dl = dl.to(u.mm)
        self.dm = dm.to(u.mm)
        self.dn = dn.to(u.mm)
    
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
        if (self.wave is not None):
            if type(self.wave.value) == np.ndarray:
                rays.addParam('Wavelength',self.wave.value)
            else:
                rays.addParam('Wavelength',np.array([self.wave.value]*len(rays)))
        if (self.order is not None):
            if type(self.order) == np.ndarray:
                rays.addParam('Order',self.order)
            else:
                rays.addParam('Order',np.array([self.order]*len(rays)))
    
    def getRays():
        '''
        Function getRays:
        Returns the Rays that have been generated by generateRays()
        '''
        if self.rays is None:
            self.generateRays()
        return self.rays

class Annulus(Source):
    
    def __init__(self,num,rin,rout,zhat=-1.,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        
        if (type(rin) != u.quantity.Quantity or type(rout) != u.quantity.Quantity):
            raise ValueError('rin and rout must be astropy Quantities')
        self.rin = rin.to(u.mm)
        self.rout = rout.to(u.mm)
        self.zhat = zhat
    
    def generateRays(self):
        rays = Rays.annulus(self.rin.value,self.rout.value,self.num,self.zhat)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays


class CircularBeam(Source):
    
    def __init__(self,num,rad,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        if (type(rad) != u.quantity.Quantity):
            raise ValueError('rad must be an astropy Quantity')
        self.rad = rad.to(u.mm)
    
    def generateRays(self):
        rays = Rays.circularbeam(self.rad.value,self.num)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays


class ConvergingBeam(Source):
    
    def __init__(self,num,zset,rin,rout,tmin,tmax,lscat,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        if (type(zset) != u.quantity.Quantity or type(rin) != u.quantity.Quantity
        or type(rout) != u.quantity.Quantity or type(tmin) != u.quantity.Quantity
        or type(tmax) != u.quantity.Quantity or type(lscat) != u.quantity.Quantity):
            raise ValueError('zset, rin, and rout must be astropy units of length. tmin, tmax, and lscat must be astropy units of angle')
        self.zset = zset.to(u.mm)
        self.rin = rin.to(t.mm)
        self.rout = rout.to(u.mm)
        self.tmin = tmin.to(u.rad)
        self.tmax = tmax.to(u.rad)
        self.lscat = lscat.to(u.arcsec)
    
    def generateRays(self):
        rays = Rays.convergingbeam(self.zset.value,self.rin.value,
        self.rout.value,self.tmin.value,self.tmax.value,self.num,
        self.lscat.value)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays


class ConvergingBeam2(Source):
    
    def __init__(self,num,zset,xmin,xmax,ymin,ymax,lscat,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        if (type(zset) != u.quantity.Quantity or type(xmin) != u.quantity.Quantity
        or type(xmax) != u.quantity.Quantity or type(ymin) != u.quantity.Quantity
        or type(ymax) != u.quantity.Quantity or type(lscat) != u.quantity.Quantity):
            raise ValueError('zset, xmin, xmax, ymin, and ymax must be astropy units of length. lscat must be an astropy unit of angle')
        self.zset = zset.to(u.mm)
        self.xmin = xmin.to(u.mm)
        self.xmax = xmax.to(u.mm)
        self.ymin = ymin.to(u.mm)
        self.ymax = ymax.to(u.mm)
        self.lscat = lscat.to(u.arcsec)
    
    def generateRays(self):
        rays = Rays.convergingbeam2(self.zset.value,self.xmin.value,
        self.xmax.value,self.ymin.value,self.ymax.value,self.num,
        self.lscat.value)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays


class PointSource(Source):
    
    def __init__(self,num,ang,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        if (type(ang) != u.quantity.Quantity):
            raise ValueError('ang must be an astropy unit of angle')
        self.ang = ang.to(u.rad)
    
    def generateRays(self):
        rays = Rays.pointsource(self.ang.value,self.num)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays


class RectArray(Source):
    
    def __init__(self,num,xsize,ysize,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        if (type(xsize) != u.quantity.Quantity or type(ysize) != u.quantity.Quantity):
            raise ValueError('xsize and ysize must be astropy units of length')
        self.xsize = xsize.to(u.mm)
        self.ysize = ysize.to(u.mm)
    
    def generateRays(self):
        rays = Rays.rectarray(self.xsize.value,self.ysize.value,self.num)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays


class RectBeam(Source):
    
    def __init__(self,num,xhalfwidth,yhalfwidth,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        if (type(xhalfwidth) != u.quantity.Quantity or 
        type(yhalfwidth) != u.quantity.Quantity):
            raise ValueError('xhalfwidth and yhalfwidth must be astropy units of length')
        self.xhalfwidth = xhalfwidth.to(u.mm)
        self.yhalfwidth = yhalfwidth.to(u.mm)
    
    def generateRays(self):
        rays = Rays.rectbeam(self.xhalfwidth.value,self.yhalfwidth.value,self.num)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays


class Subannulus(Source):
    
    def __init__(self,num,rin,rout,dphi,zhat=1.,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        if (type(rin) != u.quantity.Quantity or type(rout) != u.quantity.Quantity
        or type(dphi) != u.quantity.Quantity):
            raise ValueError('rin and rout must be astropy units of length, dphi must be an astropy unit of angle')
        self.rin = rin.to(u.mm)
        self.rout = rout.to(u.mm)
        self.dphi = dphi.to(u.rad)
        self.zhat = zhat
    
    def generateRays(self):
        rays = Rays.subannulus(self.rin.value,self.rout.value,self.dphi.value,
        self.num,self.zhat)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays


class Xslit(Source):
    
    def __init__(self,num,xin,xout,zhat=-1.,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        if (type(min) != u.quantity.Quantity or type(xout) != u.quantity.Quantity):
            raise ValueError('xin and xout must be astropy units of length')
        self.xin = xin.to(u.mm)
        self.xout = xout.to(u.mm)
        self.zhat = zhat
    
    def generateRays(self):
        rays = Rays.xslit(self.xin.value,self.xout.value,self.num,self.zhat)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays