import numpy as np
from prtp.Grating import Grating
from prtp.Combination import Combination
from prtp.Rays import Rays
import prtp.transformationsf as trans
import astropy.units as u

class GratingStack(Combination):
    '''
    Class GratingStack:
    A special kind of combination that specifically handles a group of gratings
    '''
    
    @u.quantity_input(rx=u.mm,ry=u.mm,rz=u.mm)
    def __init__(self,rx=None, ry=None, rz=None, keeporder=True):
        '''
        Initializes the GratingStack:
        
        Inputs:
        rx,ry,rz - The point about which the whole stack will rotate, see
        defineRotationPoint for more info. Must be astropy units of length
        keeporder - If True, photons will be traced to the Gratings in the order
            they were added to the stack. If False, the stack will use 
            smartTrace, where photons are sent to the nearest Grating first
        
        Notes:
        - If you want more complicated Gratings, you can modify their parameters
        using self.modifyParam(name,value), or access the Gratings themselves
        in the self.componentlist parameter
        '''
        Combination.__init__(self)
        if rx is not None:
            self.rx = rx.to(u.mm)
        if ry is not None:
            self.ry = ry.to(u.mm)
        if rz is not None:
            self.rz = rz.to(u.mm)
        self.keeporder = keeporder
    
    ## Ray-Tracing Functions
    
    def trace(self, rays, considerweights=False,eliminate='remove'):
        '''
        Function trace:
        Traces the rays through all of the Gratings in the Stack
        
        Inputs:
        rays - The rays you want to trace
        considerweights - If true, any effect that probabilistically removes
            photons will instead affect their weights
        eliminate - If 'remove', photons that miss will be removed. Otherwise,
            missed photons will be replaced with NaNs in the x-position
        
        Outputs:
        Efficiency information about the stack
        
        Notes:
        - Assumes that each Grating has been given the necessary parameters, 
        this function works with no user input.
        '''
        if self.keeporder:
            return self.defaultTrace(rays,considerweights,eliminate)
        else:
            return self.smartTrace(rays,considerweights,eliminate)
        
    
    def defaultTrace(self,rays,considerweights=False,eliminate='remove'):
        '''
        Function defaultTrace:
        Traces the Rays through the Grating Stack in the order that the Gratings
            were added to the Stack. This function will be called if 
            self.keeporder is True
        
        Inputs:
        rays - The rays you want to trace through the stack
        considerweights - Boolean saying if you want to consider the reflectivity
            of the Gratings
        eliminate - If 'remove', photons that miss will be removed. Otherwise,
            missed photons will be replaced with NaNs in the x-position
        
        Outputs:
        A tuple containing information about the efficiency of the Gratings
        '''
        # Make a blank Rays object to store the Rays that make it
        if eliminate == 'remove':
            finalrays = Rays()
        else:
            # When eliminate is NaN, we will keep track of traced rays using
            # finalrays, where np.nan represents a photon which has not yet
            # been successfully traced
            finalrays = rays.copy()
            finalrays.x[:] = np.nan
            # success stores which rays have been successfully traced.
            success = np.zeros(len(rays)).astype(bool)
        
        
        # Keep track of the input rays for when we're finished with one Grating
        inputrays = rays.copy()
        
        # Keep track of the length of the input rays
        l = rays.length(considerweights)
        
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
            if eliminate == 'remove':
                finalrays += hitrays
                # Take the rays that hit this grating out of the original Rays object
                inputrays.remove(tarray)
            else:
                # Of the rays which were just successfully traced, for which ones
                # is it their first time? Those ones need to be pulled
                rays_to_pull = np.logical_and(np.logical_not(success),tarray)
                finalrays.pullrays(rays,rays_to_pull)
                
                # Update success array
                success = np.logical_not(np.isnan(finalrays.x))
            
            # Back remaining rays up to their original position
            rays.makecopy(inputrays)
            
            if len(rays) == 0:
                break
        
        # Make it so that the original rays now contain the output
        rays.makecopy(finalrays)
        
        return ("Missed Grating Stack", l, rays.length(considerweights))
    
    
    def smartTrace(self,rays,considerweights=False,eliminate='remove'):
        '''
        Function smartTrace:
        Traces the Rays through the Grating Stack in the order that the photons
            would collide with them. This function will be called if 
            self.keeporder is False
        
        Inputs:
        rays - The rays you want to trace through the stack
        considerweights - Boolean saying if you want to consider the 
            reflectivity of the Gratings
        eliminate - If 'remove', photons that miss will be removed. Otherwise,
            missed photons will be replaced with NaNs in the x-position
        
        Outputs:
        A tuple containing information about the efficiency of the Gratings
        
        Notes:
        This function is usually slower than defaultTrace. It should only
            be used if different photons in the Rays object will encounter 
            Gratings in a different order.
        '''
        # Make a blank Rays object to store the Rays that make it
        if eliminate == 'remove':
            finalrays = Rays()
        else:
            # When eliminate is NaN, we will keep track of traced rays using
            # finalrays, where np.nan represents a photon which has not yet
            # been successfully traced
            finalrays = rays.copy()
            finalrays.x[:] = np.nan
            # success stores which rays have been successfully traced.
            success = np.zeros(len(rays)).astype(bool)
        
        # Keep track of the input rays for when we're finished with one Grating
        inputrays = rays.copy()
        
        # Keep track of the length of the input rays
        l = rays.length(considerweights)
        
        # Find the order that each photon will see the gratings
        orders = []
        for g in self.componentlist:
            orders.append(g.getDist(rays))
        # orderarr stores (for each photon) the order in which it will see
        # the gratings
        orderarr = np.stack(orders,axis=1)
        orderarr = np.argsort(orderarr).astype(float)
        
        i = 0
        while True:
            
            # Check if we've successfully traced everything
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
            
            newrays = rays.copy()
            
            self.componentlist[i].trace_to_surf(newrays)
            
            # Find which rays have hit the grating
            hit = self.componentlist[i].hit(newrays)
            
            # Keep those which have hit the Grating
            if eliminate == 'remove':
                newrays.remove(np.logical_or(np.logical_not(hit),np.logical_not(tarray)))
            else:
                newrays.x[np.logical_or(np.logical_not(hit),np.logical_not(tarray))] = np.nan
            
            if (np.sum(hit) != 0):
                # These operations can only be done on non-empty Rays Objects
                
                # Trace and save the Rays which hit the Grating
                self.componentlist[i].trace(newrays,eliminate=eliminate)
                
                # Add the hitrays to our final tally
                if eliminate == 'remove':
                    finalrays += newrays
                    # Take the rays that hit this grating out of the original Rays object
                    inputrays.remove(hit)
                else:
                    # Of the rays which were just successfully traced, for which ones
                    # is it their first time? Those ones need to be pulled
                    rays_to_pull = np.logical_and(np.logical_and(np.logical_not(success),hit),tarray)
                    finalrays.pullrays(newrays,rays_to_pull)
                    
                    # Update success array
                    success = np.logical_not(np.isnan(finalrays.x))
                
            
            # Use this hit trutharray to find which of the original rays have 
            # hit the Grating
            test = tarray.copy()
            tarray[tarray] = hit[tarray]
            
            # Remove the rays which have hit this Grating
            if eliminate == 'remove':
                rays.remove(tarray)
            else:
                rays.x[tarray] = np.nan
            
            # Update the orderarr to rotate out those which needed to be traced
            # to this Grating
            orderarr[test] = np.roll(orderarr[test],-1,axis=1)
            
            # Set the already tried indices to -1 to make sure we don't try them
            # again
            orderarr[test,-1] = -1
            
            # Update orderarr to remove the hit photons
            if eliminate == 'remove':
                orderarr = np.delete(orderarr,np.where(tarray)[0],0)
            else:
                orderarr[tarray] == np.nan
            
            # Go to the next Grating
            i += 1
            
            if (i >= len(self.componentlist)):
                i = 0
        
        # Make it so that the original rays now contain the output
        rays.makecopy(finalrays)
        
        return ("Missed Grating Stack", l, rays.length(considerweights))