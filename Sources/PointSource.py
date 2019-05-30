import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class PointSource(Source):
    
    def __init__(self,num,ang,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.ang = ang
    
    def generateRays(self):
        rays = Rays.pointsource(self.ang,self.num)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays