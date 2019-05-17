import numpy as np
import matplotlib.pyplot as plt
from prtp.Rays import Rays
import prtp.specialFunctions as func
import prtp.transformationsf as trans
from prtp.WolterOptic import WolterOptic

class WolterTypeOne(WolterOptic):
    
    def __init__(self,x=0,y=0,z=0,nx=0,ny=0,nz=1,r0=1,z0=1,psi=1):
        '''
        WolterTypeOne Object:
        - Combines a WolterPrimary and a WolterSecondary Object into one single object
        The WolterTypeOne Object will process rays faster than you would with the two individual objects, owing to the fact that the rays only have to be transformed once.
        DO NOT use a WolterTypeOne object if:
            - You need to move the two components relative to each other so that they no longer share a common focus (e.g: you try to handle misalignments in the two optics separately)
            - You need to perform a tracing more complicated than:
                - Trace to Primary
                - Reflect
                - Trace to Secondary
            e.g: with WolterTypeOne Objects you cannot add scattering on the surface of the Primary (since the rays are automatically traced to the secondary, giving the user no time to perform further modifications)
        
        
        Parameters:
        x,y,z - The Cartesian Coordinates of the focus of the optic
        nx,ny,nz - The components of the vector pointing outwards from the focus of the optic
        r0 - The radius of the optic when it converges with the WolterSecondary
        z0 - The position along the n-direction at which this optic converges with the WolterSecondary (the focus is at position 0)
        psi - Some ratio of angles, I don't know, I think its always 1
        '''
        WolterOptic.__init__(self,x,y,z,nx,ny,nz,r0,z0,psi)


    ## Tracing Rays to the Optic:
    
    def tracefunction(self,rays,autoreflect):
        '''
        Function tracefunction:
        Determines what happens to the rays once they has been transformed, should only be called by the WolterOptic superclass's trace function
        '''
        rays.wolterprimary(self.r0,self.z0,self.psi)
        rays.reflect()
        rays.woltersecondary(self.r0,self.z0,self.psi)
        
        if autoreflect:
            rays.reflect()
    
    def trace(self,rays,autoreflect=True,remove=False):
        '''
        Function trace:
        Traces rays to the wolter optic
        
        Inputs:
        rays - The Rays object you want to trace to the optic
        autoreflect - If True, the photons will be automatically reflected off of the secondary. If False, the photons will not be reflected off the secondary optic
        remove - If True, photons that missed the optic will be automatically removed. If False, photons that missed will be given nans in every quantity
        
        Outputs:
        A tuple containing the original number of rays and the final number of rays (that hit the optic)
        '''
        return self.tracehelper(rays,self.tracefunction,autoreflect,remove)