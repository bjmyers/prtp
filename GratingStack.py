import numpy as np
from prtp.Grating import Grating
from prtp.Combination import Combination
from prtp.Rays import Rays
import prtp.transformationsf as trans

class GratingStack(Combination):
    '''
    Class GratingStack:
    A special kind of combination that specifically handles a group of gratings
    '''
    
    def __init__(self,radial=True,d=160,order=0,wave=100,autoreflect=False,
    rx=0, ry=0, rz=0):
        '''
        Initializes the GratingStack:
        
        Inputs:
        radial, if True, all Gratings in this stack will call their radgrat 
        function if False, all Gratings in this stack will call their grat function
        d - Grating period, can be an integer of an array the same length as the 
        incoming Rays
        order - The order of the incoming photons, can be an integer of an array 
        the same length as the incoming Rays
        wave - The wavelength of the incoming photons, can be an integer of an 
        array the same length as the incoming Rays
        autoreflect - If True, Rays will automatically be reflected off of the
        Grating they hit
        rx,ry,rz - The point about which the whole stack will rotate, see
        defineRotationPoint for more info
        
        Notes:
        - Currently, GratingStack does not support gratings of different types
        or different periods being grouped together
        '''
        Combination.__init__(self)
        self.radial = radial
        self.d = d
        self.order = order
        self.wave = wave
        self.autoreflect = True
        self.rx = rx
        self.ry = ry
        self.rz = rz
    
    
    def defineRotationPoint(self,x=0,y=0,z=0):
        '''
        Function defineRotationPoint:
        Defines the point about which the entire Stack can rotate
        
        Inputs:
        x,y,z - The Coordinates of the Rotation Point
        
        Outputs:
        None
        
        Notes: 
        Calling this function is required if you try to rotate the whole stack.
        Using Combination.applyToAll and a rotation function will rotate each
        Grating about their centers, rather than rotating the whole stack about
        this point
        '''
        self.rx = x
        self.ry = y
        self.rz = z
    
    
    def unitrotate(self,theta=0,axis=1):
        '''
        Function unitrotate:
        Rotates the entire Graint Stack about the rotationpoint and about a
        unit axis
        
        Inputs:
        theta - The amount by which you want to rotate
        axis - integer input of 1, 2, or 3 to rotate about the x, y, or z axes, respectively.
        
        Outputs:
        None
        '''
        
        for g in self.componentlist:
            
            # Rotates the Grating's two vectors
            g.unitrotate(theta,axis)
            
            # Move the Grating's so the rotation point is about the origin 
            g.translate(-self.rx,-self.ry,-self.rz)
            
            # Rotate the Grating's position
            g.x,g.y,g.z = trans.rotatevector(g.x,g.y,g.z,theta,axis)
            
            # Translate the origin back down
            g.translate(self.rx,self.ry,self.rz)
    
    
    def rotate(self,theta,ux,uy,uz):
        '''
        Function unitrotate:
        Rotates the entire Graint Stack about the rotationpoint and about a
        user-defined axis
        
        Inputs:
        theta - The amount by which you want to rotate
        ux,uy,uz - The x, y, and z components of the vector about which you want 
        to rotate
        
        Outputs:
        None
        '''
        
        for g in self.componentlist:
            
            # Rotates the Grating's two vectors
            g.rotate(theta,ux,uy,uz)
            
            # Move the Grating's so the rotation point is about the origin 
            g.translate(-self.rx,-self.ry,-self.rz)
            
            # Rotate the Grating's position
            g.x,g.y,g.z,q1,q2,q3 = trans.rotateaxis(g.x,g.y,g.z,ux,uy,uz,theta)
            
            # Translate the origin back down
            g.translate(self.rx,self.ry,self.rz)
    
    def translate(self,dx=0,dy=0,dz=0):
        '''
        Function translate
        Translates the GratingStack in three-dimensions
        
        Inputs:
        dx,dy,dz - The amount to move in x, y, and z, respectively
        
        Outputs:
        None
        
        Notes:
        - This move is relative, not absolute. That is, you will move BY dx, dy, and z, you will not move TO dx, dy, and dz
        '''
        self.applyToAll(Grating.translate,dx=dx,dy=dy,dz=dz)
    
    
    def trace(self, rays):
        
        # Make a blank Rays object to store the Rays that make it
        finalrays = Rays()
        
        # Keep track of the input rays for when we're finished with one Grating
        inputrays = rays
        
        # Keep track of how many photons go onto each grating
        effs = []
        
        # Iterate through each Grating Object
        for g in self.componentlist:
            # Through each pass we need to ensure that the rays that make it are 
            # placed into a final rays object
            # All those that miss are passed to the next Grating
            g.trace_to_surf(rays)
            
            # Find which rays have hit the grating
            tarray = g.hit(rays)
            
            hitrays = rays.split(tarray)
            
            if self.radial:
                g.radgrat(hitrays,self.d,self.order,self.wave,self.autoreflect)
            else:
                g.grat(hitrays,self.d,self.order,self.wave,self.autoreflect)
            
            # Add the hitrays to our final tally
            finalrays += hitrays
            
            # Record how many photons hit this Grating Object
            effs.append(len(hitrays))
            
            # Take the rays that hit this grating out of the original Rays object
            inputrays.remove(tarray)
            
            # Back remaining rays up to their original position
            rays = inputrays
            
            if len(rays) == 0:
                break
        
        self.effs = effs
        
        return finalrays
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        