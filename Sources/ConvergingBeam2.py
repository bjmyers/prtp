import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

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