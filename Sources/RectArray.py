import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class RectArray(Source):
    
    def __init__(self,num,xsize,ysize,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.xsize = xsize
        self.ysize = ysize
    
    def generateRays(self):
        rays = Rays.rectarray(self.xsize,self.ysize,self.num)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays