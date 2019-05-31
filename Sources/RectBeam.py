import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source
import astropy.units as u

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