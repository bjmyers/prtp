import numpy as np
from prtp.Rays import Rays
import astropy.units as u
import prtp.transformationsf as trans

class Combination:
    '''
    Class Combination:
    A combination object is a group of several components that the user wants to group together
    When Rays are traced to a Combination Object, they will be traced to each Component individually and only those components who missed every component will be removed
    Tracing to Combinations will also give detailed information about how each Component affected the overall efficiency
    '''
    
    def __init__(self):
        self.componentlist = []
    
    
    def addComponent(self,comp,index=None):
        '''
        Function addComponent:
        Adds a Component to the componentlist
        
        Inputs:
        comp - The component to add
        index - Where in the list you want to add the component, if None, it will be added to the end
        '''
        if index is None:
            self.componentlist.append(comp)
        else:
            self.componentlist.insert(index,comp)
    
    def applyToAll(self, func, **kwargs):
        '''
        Function applyToAll:
        Applies a function to each component in the Combination Object
        (e.g: applyToAll(self,FlatComponent.UnitRotate,theta=.1,axis=2) will rotate each component in the Combination .1 radians about the y-axis
        
        Inputs:
        func - The function you want to apply to each element
        **kwargs - Any other arguments that func requires
        
        Outputs:
        None
        '''
        for c in self.componentlist:
            func(c,**kwargs)
    
    ## Movement Functions:
    
    @u.quantity_input(x=u.mm,y=u.mm,z=u.mm)
    def defineRotationPoint(self,x=0*u.mm,y=0*u.mm,z=0*u.mm):
        '''
        Function defineRotationPoint:
        Defines the point about which the entire Stack can rotate
        
        Inputs:
        x,y,z - The Coordinates of the Rotation Point, must be astropy units of
            length
        
        Outputs:
        None
        
        Notes: 
        Calling this function is required if you try to rotate the whole stack.
        Using Combination.applyToAll and a rotation function will rotate each
        Grating about their centers, rather than rotating the whole stack about
        this point
        '''
        self.rx = x.to(u.mm)
        self.ry = y.to(u.mm)
        self.rz = z.to(u.mm)
    
    @u.quantity_input(theta=u.rad)
    def unitrotate(self,theta=0*u.rad,axis=1):
        '''
        Function unitrotate:
        Rotates the entire Graint Stack about the rotationpoint and about a
        unit axis
        
        Inputs:
        theta - The amount by which you want to rotate, must be an astropy unit
            of angle
        axis - integer input of 1, 2, or 3 to rotate about the x, y, or z axes, respectively.
        
        Outputs:
        None
        '''
        theta = theta.to(u.rad)
        
        for g in self.componentlist:
            
            # Rotates the Grating's two vectors
            g.unitrotate(theta,axis)
            
            # Move the Grating's so the rotation point is about the origin 
            g.translate(-self.rx,-self.ry,-self.rz)
            
            # Rotate the Grating's position
            g.x,g.y,g.z = trans.rotatevector(g.x,g.y,g.z,theta,axis)
            
            # Translate the origin back down
            g.translate(self.rx,self.ry,self.rz)
    
    @u.quantity_input(theta=u.rad)
    def rotate(self,theta=0*u.rad,ux=1,uy=0,uz=0):
        '''
        Function unitrotate:
        Rotates the entire Graint Stack about the rotationpoint and about a
        user-defined axis
        
        Inputs:
        theta - The amount by which you want to rotate, must be an astropy unit
            of angle
        ux,uy,uz - The x, y, and z components of the vector about which you want 
        to rotate
        
        Outputs:
        None
        '''
        theta = theta.to(u.rad)
        
        for g in self.componentlist:
            
            # Rotates the Grating's two vectors
            g.rotate(theta,ux,uy,uz)
            
            # Move the Grating's so the rotation point is about the origin 
            g.translate(-self.rx,-self.ry,-self.rz)
            
            # Rotate the Grating's position
            g.x,g.y,g.z,q1,q2,q3 = trans.rotateaxis(g.x.value,g.y.value,g.z.value,ux,uy,uz,theta.value)
            # Restore units
            g.x *= u.mm
            g.y *= u.mm
            g.z *= u.mm
            
            # Translate the origin back down
            g.translate(self.rx,self.ry,self.rz)
    
    @u.quantity_input(x=u.mm,y=u.mm,z=u.mm)
    def translate(self,dx=0*u.mm,dy=0*u.mm,dz=0*u.mm):
        '''
        Function translate
        Translates the GratingStack in three-dimensions
        
        Inputs:
        dx,dy,dz - The amount to move in x, y, and z, respectively, must be
            astropy units of length
        
        Outputs:
        None
        
        Notes:
        - This move is relative, not absolute. That is, you will move BY dx, dy, and z, you will not move TO dx, dy, and dz
        '''
        for g in self.componentlist:
            g.translate(dx,dy,dz)