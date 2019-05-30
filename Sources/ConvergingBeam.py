import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class ConvergingBeam(Source):
    
    def __init__(self,num,zset,rin,rout,tmin,tmax,lscat,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.zset = zset
        self.rin = rin
        self.rout = rout
        self.tmin = tmin
        self.tmax = tmax
        self.lscat = lscat
    
    def generateRays(self):
        rays = Rays.convergingbeam(self.zset,self.rin,self.rout,self.tmin,self.tmax,self.num,self.lscat)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays