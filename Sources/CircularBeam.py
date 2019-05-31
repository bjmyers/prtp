import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source
import astropy.units as u

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