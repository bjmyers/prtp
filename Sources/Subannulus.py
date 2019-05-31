import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source
import astropy.units as u

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