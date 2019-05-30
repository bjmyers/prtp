import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class Annulus(Source):
    
    def __init__(self,num,rin,rout,zhat=-1.,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.rin = rin
        self.rout = rout
        self.zhat = zhat
    
    def generateRays(self):
        rays = Rays.annulus(self.rin,self.rout,self.num,self.zhat)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays