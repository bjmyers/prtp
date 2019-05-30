import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class Xslit(Source):
    
    def __init__(self,num,xin,xout,zhat=-1.,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.xin = xin
        self.xout = xout
        self.zhat = zhat
    
    def generateRays(self):
        rays = Rays.xslit(self.xin,self.xout,self.num,self.zhat)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays