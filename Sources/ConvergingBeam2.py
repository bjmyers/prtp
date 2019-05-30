import numpy as np
from prtp.Rays import Rays
from prtp.Sources.Source import Source

class ConvergingBeam2(Source):
    
    def __init__(self,num,zset,rin,rout,tmin,tmax,lscat,wave=None,order=None):
        Source.__init__(self,num,wave,order)
        self.zset = zset
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.lscat = lscat
    
    def generateRays(self):
        rays = Rays.convergingbeam2(self.zset,self.xmin,self.xmax,self.ymin,self.ymax,self.num,self.lscat)
        self.addParams(rays)
        self.misalign(rays)
        self.rays = rays
        return rays