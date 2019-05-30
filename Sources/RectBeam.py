import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class ConvergingBeam(Source):
    
    def __init__(self,num,xhaflwidth,yhalfwidth,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.xhalfwidth = xhalfwidth
        self.yhalfwidth = yhalfwidth
    
    def generateRays(self):
        rays = Rays.rectbeam(self.xhalfwidth,self.yhalfwidth,self.num)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays