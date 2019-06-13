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
        axis - integer input of 1, 2, or 3 to rotate about the x, y, or z axes, 
            respectively.
        
        Outputs:
        None
        
        Notes:
        - Theta can be a single value or an array/list. If it is array-like, it
            must have the same length as the component list as it specifies the
            rotation angle of each component
        '''
        theta = theta.to(u.rad)
        
        try:
            len(theta)
            if (len(theta) != len(self.componentlist)):
                raise ValueError('Theta Array must be the same length as the component list')
            else:
                arr = True
        except:
            arr = False
        # arr is True if theta is an array, False otherwise
        
        if (arr):
            for i in range(len(self.componentlist)):
                g = self.componentlist[i]
                # Rotates the Grating's two vectors
                g.unitrotate(theta[i],axis)
                # Move the Grating's so the rotation point is about the origin 
                g.translate(-self.rx,-self.ry,-self.rz)
                # Rotate the Grating's position
                g.x,g.y,g.z = trans.rotatevector(g.x,g.y,g.z,theta[i],axis)
                # Translate the origin back down
                g.translate(self.rx,self.ry,self.rz)
        else:
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
        
        Notes:
        - Theta can be a single value or an array/list. If it is array-like, it
            must have the same length as the component list as it specifies the
            rotation angle of each component
        '''
        theta = theta.to(u.rad)
        
        try:
            len(theta)
            if (len(theta) != len(self.componentlist)):
                raise ValueError('Theta Array must be the same length as the component list')
            else:
                arr = True
        except:
            arr = False
        # arr is True if theta is an array, False otherwise
        
        if (arr):
            for i in range(len(self.componentlist)):
                g = self.componentlist[i]
                # Rotates the Grating's two vectors
                g.rotate(theta[i],ux,uy,uz)
                # Move the Grating's so the rotation point is about the origin 
                g.translate(-self.rx,-self.ry,-self.rz)
                # Rotate the Grating's position
                g.x,g.y,g.z,q1,q2,q3 = trans.rotateaxis(g.x.value,g.y.value,g.z.value,ux,uy,uz,theta[i].value)
                # Restore units
                g.x *= u.mm
                g.y *= u.mm
                g.z *= u.mm
                # Translate the origin back down
                g.translate(self.rx,self.ry,self.rz)
        else:
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
    
    
    @u.quantity_input(theta=u.rad)
    def unitrotateinplace(self,theta,axis=1):
        '''
        Function unitrotateinplace:
        Rotates each of the components in this combination about a unit axis.
        Unlike unitrotate which rotates each component about a single point,
        this function rotates each component about their center. Leaving them
        "in place" and not changing the position of their centers.
        
        Inputs:
        theta - The amount by which you want to rotate, must be an astropy unit
            of angle
        axis - integer input of 1, 2, or 3 to rotate about the x, y, or z axes, 
            respectively.
        
        Outputs:
        None
        
        Notes:
        - Theta must be an array of quantities the same length as self.componentlist
            If you wished to rotate each component by the same amount, use
            applyToAll()
        '''
        if (len(theta) != len(self.componentlist)):
            raise ValueError('Theta Array must be the same length as the component list')
        for i in range(len(self.componentlist)):
            self.componentlist[i].unitrotate(theta[i],axis)
    
    
    @u.quantity_input(theta=u.rad)
    def rotateinplace(self,theta,ux=1,uy=0,uz=0):
        '''
        Function unitrotateinplace:
        Rotates each of the components in this combination about a user-defined
        axis. Unlike unitrotate which rotates each component about a single point,
        this function rotates each component about their center. Leaving them
        "in place" and not changing the position of their centers.
        
        Inputs:
        theta - The amount by which you want to rotate, must be an astropy unit
            of angle
        ux,uy,uz - Components of the vector about which you want to rotate, does
            not need to be a unit vector (can have any length)
        
        Outputs:
        None
        
        Notes:
        - Theta must be an array of quantities the same length as self.componentlist
            If you wished to rotate each component by the same amount, use
            applyToAll()
        '''
        if (len(theta) != len(self.componentlist)):
            raise ValueError('Theta Array must be the same length as the component list')
        for i in range(len(self.componentlist)):
            self.componentlist[i].rotate(theta[i],ux,uy,uz)
    
    
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
    
    
    ## Trace Function:
    
    def trace(self,rays,considerweights=False,eliminate='remove'):
        '''
        Function trace:
        Traces rays to each component in the Combination, keeps the photons 
        which hit either component
        
        Inputs:
        rays - The Rays object you want to trace
        considerweights - If True, the weighting of rays will be used when
            trying to remove rays based on a probability, with reflectivity,
            for example.
        eliminate - If 'remove', photons that miss will be removed from the 
            Rays object. If it is anything else, the photons which miss will
            be given NaN for their x-position
        
        Outputs:
        Efficiency information about the combination
        '''
        # Store the length of the incoming rays
        l = rays.length(considerweights)
        
        # Initialize a list in which efficiencies will be stored
        effs = []
        
        # Make a blank Rays object to store the Rays that make it
        finalrays = Rays()
        
        # Keep track of the input rays for when we're finished with one Mirror
        inputrays = rays.copy()
        temprays = rays.copy()
        
        # Iterate through each Mirror Object
        for r in self.componentlist:
            # Through each pass we need to ensure that the rays that make it are 
            # placed into a final rays object
            # All those that miss are passed to the next Mirror
            eff = r.trace(temprays,considerweights=considerweights,eliminate='nan')
            
            # Some components return a list of efficiencies, check if this is the case
            if type(eff) == list:
                effs.extend(eff)
            else:
                effs.append(eff)
            
            # Find the Rays which missed the 
            tarray = np.logical_not(np.isnan(temprays.x))
            hitrays = temprays.split(tarray)
            
            # Take the rays that hit this grating out of the original Rays object
            inputrays.remove(tarray)
            
            # Back remaining rays up to their original position
            temprays = inputrays.copy()

            # Make sure at least some rays have hit the mirror
            if (len(hitrays) == 0):
                continue
            
            # Add the hitrays to our final tally
            finalrays += hitrays
            
            # If there are no rays left, we can stop
            if len(temprays) == 0:
                break
        
        # Make it so that the original rays now contain the output
        rays.makecopy(finalrays)
        
        # Sort efficiency outputs into rays lost by missing the components and
        # rays lost through other effects

        # hit_but_lost stores the number of rays which hit a component but were
        # lost through other effects
        hit_but_lost = 0
        for e in effs:
            # All missed efficiencies start with "Missed"
            if e[0][:6] != 'Missed':
                hit_but_lost += e[2] - e[1]
        
        eff1 = ('Missed Combination',l,rays.length(considerweights)+hit_but_lost)
        
        eff2 = ('Lost By Combination - Other',
        rays.length(considerweights)+hit_but_lost,rays.length(considerweights))
        
        return [eff1,eff2]
