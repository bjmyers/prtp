import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
from prtp.FlatComponent import FlatComponent
import prtp.transformationsf as trans

class Grating(FlatComponent):
    
    ## Initialization Function
    
    def __init__(self,x=0,y=0,z=0,nx=0,ny=0,nz=1,fx=0,fy=1,fz=0,l=1,w=1):
        '''
        Initializes a Grating Object, requires the following arguments:
        
        x,y,z - The position of a point along the grating
        nx,ny,nz - The components of a vector normal to the grating surface
        fx,fy,fz - The components of a vector pointing towards the grating focus
        l - The length of the grating, that is, its extent in the f direction
        w - The width of the grating, that is, its extent in the fxn direction
        '''
        FlatComponent.__init__(self,x,y,z,nx,ny,nz,fx,fy,fz)
        self.l = l
        self.w = w
    
    ## Analyzing Traced Rays
    
    def hit(self,rays):
        '''
        Function hit:
        Given Rays that have been traced to the Component, returns a trutharray detailing which photons have fallen within the rectangular region of this grating.
        
        Inputs:
        rays - A Rays object that has been traced to the Component Plane
        
        Outputs:
        tarray - A trutharray, containing True if the photon has fallen within the grating's dimensions or False if the photon missed the grating
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        - Make sure the length and the width have been defined when the Grating was defined
        '''
        x,y = self.getPosns(rays)
        return np.logical_and(np.abs(x) < self.w/2, np.abs(y) < self.l/2)
    
    def removemissed(self,rays):
        '''
        Function removemissed:
        Given Rays that have been traced to the Component, returns the photons which have hit the Grating.
        
        Inputs:
        rays - A Rays object that has been traced to the Component Plane
        
        Outputs:
        rays - The updated Rays object, containing only those that have hit the Grating
        
        Notes:
        - The function cannot tell if the Rays have been traced, so this is up to the user.
        - Make sure the length and the width have been defined when the Grating was defined
        '''
        tarray = self.hit(rays)
        rays.remove(tarray)
    
    ## Grating Functions:
    # Handle reflecting off of the gratings
    
    def grat(self, rays, d=160, order=0, wave=100):
        vel = [rays.l,rays.m,rays.n]
        # Define normal vector, focus vector, and their cross product, stacked on themselves so they can be dotted into the rays
        norm = np.tile([self.nx,self.ny,self.nz],(len(rays),1)).transpose()
        focus = np.tile([self.sx,self.sy,self.sz],(len(rays),1)).transpose()
        sxn = np.cross([self.sx,self.sy,self.sz],[self.nx,self.ny,self.nz])
        sxn = np.tile(sxn,(len(rays),1)).transpose()
        
        # Define velocity components aligned with Grating (if it were in the xy-plane) using dot products
        gratl = (vel * sxn).sum(0)
        gratm = (vel * focus).sum(0)
        gratn = (vel * norm).sum(0)
        
        l,m,n = trans.grat(gratl,gratm,gratn,d,order,wave)
        rays.set(l=l,m=m,n=n)
    
    
    
    
    
    
    #TODO: Make groove density defining functions, figure out how to get rays to reflect off gratings

## Testing Code:
# r = Rays.pointsource(5,1000)
# g = Grating(0,0,5,0,0,1,0,1,0,5,5)
# g.rotate(.01,1,1,1)
# g.translate(0,0,1)
# r = g.trace_to_surf(r)
# x,y = g.getPosns(r)
# 
# # r.scatter3d()
# # plt.figure()
# # plt.scatter(x,y)
# # plt.show()
# 
# # a = g.hit(r)
# # g.removemissed(r)
# # r.scatter3d()
# 
# g.grat(r)

















