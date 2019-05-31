import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source
import astropy.units as u

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