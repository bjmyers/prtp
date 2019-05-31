import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source
import astropy.units as u

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