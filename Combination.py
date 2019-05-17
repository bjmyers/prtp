import numpy as np
from prtp.Rays import Rays

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