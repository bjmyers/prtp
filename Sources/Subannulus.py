import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class Subannulus(Source):
    
    def __init__(self,num,rin,rout,dphi,zhat=1.,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.rin = rin
        self.rout = rout
        self.dphi = dphi
        self.zhat = zhat
    
    def generateRays(self):
        rays = Rays.subannulus(self.rin,self.rout,self.dphi,self.num,self.zhat)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays