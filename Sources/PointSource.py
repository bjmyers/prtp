import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source
import astropy.units as u

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