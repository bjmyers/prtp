import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source
import astropy.units as u

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