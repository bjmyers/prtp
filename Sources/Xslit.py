import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source
import astropy.units as u

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