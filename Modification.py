import numpy as np

class Modification:
    
    def __init__(self,function=None):
        self.function = function
    
    
    def trace(self,rays,considerweights=False):
        self.function(rays,considerweights)
        return