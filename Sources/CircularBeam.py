import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class CircularBeam(Source):
    
    def __init__(self,num,rad,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.rad = rad
    
    def generateRays(self):
        rays = Rays.circularbeam(self.rad,self.num)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays