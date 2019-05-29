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
    rx=0, ry=0, rz=0, keeporder=True):
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
        keeporder - If True, photons will be traced to the Gratings in the order
            they were added to the stack. If False, the stack will use 
            smartTrace, where photons are sent to the nearest Grating first
        
        Notes:
        - All parameters except rx,ry,rz are used in the most basic form of
        GratingStack, one where all Gratings have the same properties. If this
        is the type of GratingStack you want, make sure all Gratings have these
        parameters by calling self.setDefaultParams()
        - If you want more complicated Gratings, you can modify their parameters
        using self.modifyParam(name,value), or access the Gratings themselves
        in the self.componentlist parameter
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
        self.keeporder = keeporder
    
    
    ## Parameter Functions
    # The parameters of a Grating are very important for its behavior, these
    # functions allow the user to modify the parameters of every Grating in the
    # stack
    
    def setDefaultParams():
        '''
        Function setDefaultParams:
        Sets the parameters of every Grating in this stack to the ones defined
        when the GratingStack was initialized
        
        Inputs:
        None
        
        Outputs:
        None
        '''
        for g in self.componentlist:
            g.d = self.d
            g.order = self.order
            g.wave = self.order
    
    
    def modifyParam(name,value):
        '''
        Function modifyParam:
        Changes a parameter of each Grating in the Stack
        
        Inputs:
        name - String, the name of the parameter you want to add or change
        value - The value you want to assign to that parameter
        
        Notes:
        - This function will assign the same value to each Grating in the Stack
        '''
        for g in self.componentlist:
            setattr(g,name,value)
    
    
    ## Movement Functions
    
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
    
    
    ## Ray-Tracing Functions
    
    def trace(self, rays, considerweights=False):
        '''
        Function trace:
        Traces the rays through all of the Gratings in the Stack
        
        Inputs:
        rays - The rays you want to trace
        
        Outputs:
        finalrays - The rays that have been traced through the Stack
        
        Notes:
        - Assumes that each Grating has been given the necessary parameters, 
        this function works with no user input.
        '''
        if self.keeporder:
            return self.defaultTrace(rays,considerweights)
        else:
            return self.smartTrace(rays,considerweights)
        
        
    
    def defaultTrace(self,rays,considerweights=False):
        '''
        Function defaultTrace:
        Traces the Rays through the Grating Stack in the order that the Gratings
            were added to the Stack. This function will be called if 
            self.keeporder is True
        
        Inputs:
        rays - The rays you want to trace through the stack
        considerweights - Boolean saying if you want to consider the reflectivity
            of the Gratings
        
        Outputs:
        A tuple containing information about the efficiency of the Gratings
        '''
        # Make a blank Rays object to store the Rays that make it
        finalrays = Rays()
        
        # Keep track of the input rays for when we're finished with one Grating
        inputrays = rays
        
        # Keep track of the length of the input rays
        l = len(rays)
        
        # Iterate through each Grating Object
        for g in self.componentlist:
            # Through each pass we need to ensure that the rays that make it are 
            # placed into a final rays object
            # All those that miss are passed to the next Grating
            g.trace_to_surf(rays)
            
            # Find which rays have hit the grating
            tarray = g.hit(rays)
            
            hitrays = rays.split(tarray)
            
            # Make sure at least some rays have hit the grating
            if (len(hitrays) == 0):
                continue
            
            g.trace(hitrays,considerweights)
            
            # Add the hitrays to our final tally
            finalrays += hitrays
            
            # Take the rays that hit this grating out of the original Rays object
            inputrays.remove(tarray)
            
            # Back remaining rays up to their original position
            rays = inputrays
            
            if len(rays) == 0:
                break
        
        # Make it so that the original rays now contain the output
        rays.makecopy(finalrays)
        
        return ("Missed Grating Stack", l, len(rays))
    
    
    def smartTrace(self,rays,considerweights=False):
        '''
        Function smartTrace:
        Traces the Rays through the Grating Stack in the order that the photons
            would collide with them. This function will be called if 
            self.keeporder is False
        
        Inputs:
        rays - The rays you want to trace through the stack
        considerweights - Boolean saying if you want to consider the 
            reflectivity of the Gratings
        
        Outputs:
        A tuple containing information about the efficiency of the Gratings
        
        Notes:
        This function is usually slower than defaultTrace. It should only
            be used if different photons in the Rays object will encounter 
            Gratings in a different order.
        '''
        # Make a blank Rays object to store the Rays that make it
        finalrays = Rays()
        
        # Keep track of the input rays for when we're finished with one Grating
        inputrays = rays
        
        # Keep track of the length of the input rays
        l = len(rays)
        
        # Find the order that each photon will see the gratings
        orders = []
        for g in self.componentlist:
            orders.append(g.getDist(rays))
        # orderarr stores (for each photon) the order in which it will see
        # the gratings
        orderarr = np.stack(orders,axis=1)
        orderarr = np.argsort(orderarr)
        
        i = 0
        while True:
            
            # Check if we've gotten everything
            if (orderarr[:,0] == -1).all():
                break
            
            # Find which rays need to be trace to this Grating
            tarray = (orderarr[:,0] == i)
            
            if (np.sum(tarray) == 0):
                # Go to the next Grating
                i += 1
                if (i >= len(self.componentlist)):
                    i = 0
                
                continue 
            
            newrays = rays.split(tarray)
            
            self.componentlist[i].trace_to_surf(newrays)
            
            # Find which rays have hit the grating
            hit = self.componentlist[i].hit(newrays)
            
            # Keep those which have hit the Grating
            newrays.remove(np.logical_not(hit))
            
            if (np.sum(hit) != 0):
                # These operations can only be done on non-empty Rays Objects
                
                # Trace and save the Rays which hit the Grating
                self.componentlist[i].trace(newrays)
                finalrays += newrays
            
            # Use this hit trutharray to find which of the original rays have 
            # hit the Grating
            test = tarray.copy()
            tarray[tarray] = hit
            
            # Remove the rays which have hit this Grating
            rays.remove(tarray)
            
            # Update the orderarr to rotate out those which needed to be traced
            # to this Grating
            orderarr[test] = np.roll(orderarr[test],-1,axis=1)
            
            # Set the already tried indices to -1 to make sure we don't try them
            # again
            orderarr[test,-1] = -1
            
            # Update orderarr to remove the hit photons
            orderarr = np.delete(orderarr,np.where(tarray)[0],0)
            
            # Go to the next Grating
            i += 1
            
            if (i >= len(self.componentlist)):
                i = 0

        
        # Make it so that the original rays now contain the output
        rays.makecopy(finalrays)
        
        return ("Missed Grating Stack", l, len(rays))