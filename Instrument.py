import numpy as np


class Instrument:
    # The main class for simulating optical systems, each Instrument object
    # will contain a list of components you want to trace the rays to in order
    
    def __init__(self,source = None, waves = 1., orders = 0, weighting = False):
        '''
        Initialization Function:
        Creates the Instrument Objest
        
        Inputs:
        source - A Rays Object, these are the photons which will be sent through
            the Instrument
        waves - The wavelength of the photons. Can be an int/float if you wish
            for every photon to have the same wavelength. Can also be an array
            the same length as the Rays object if you wish to specify the 
            wavelength of each photon
        orders - The order of the photons. Can be an int/float if you wish
            for every photon to have the same order. Can also be an array
            the same length as the Rays object if you wish to specify the 
            order of each photon
        
        Notes: 
        - componentlist is the list of all the components you will trace 
            photons to, in order.
        - effs is the list of efficiencies for each component, it will be filled
            when photons are traced to each component
        '''
        self.source = source
        
        # Make sure that waves and sources are arrays with the same length 
        # as source
        if not (type(waves) == list or type(waves) == np.ndarray):
            waves = np.array([waves] * len(source))
        if not (type(orders) == list or type(orders) == np.ndarray):
            orders = np.array([orders] * len(source))
            
        # Add waves and orders to the Rays object as parameters
        # If the parameters already exist, replace the values
        if self.source.getParam('Wavelength') is None:
            self.source.addParam('Wavelength',waves)
        else:
            self.source.removeParam('Wavelength')
            self.source.addParam('Wavelength',waves)
        if self.source.getParam('Order') is None:
            self.source.addParam('Order',orders)
        else:
            self.source.removeParam('Order')
            self.source.addParam('Order',orders)
        
        self.weighting = weighting
        if weighting:
            self.source.addWeighting()
        
        self.componentlist = []
        self.effs = []
    
    
    ## Changing Components
    
    def addComponent(self,comp,index=None):
        '''
        Function addComponent:
        Adds a Component to the componentlist
        
        Inputs:
        comp - The component to add
        index - Where in the list you want to add the component, if None, it 
            will be added to the end
        '''
        if index is None:
            self.componentlist.append(comp)
        else:
            self.componentlist.insert(index,comp)
    
    def removeComponent(self,index):
        '''
        Function removeComponent:
        Removes a component from the component list
        
        Inputs:
        index - The index of the component you want to remove
        
        Outputs:
        The component that was removed
        '''
        return self.componentlist.pop(index)
    
    
    ## Simulation Functions
    
    def simulate(self):
        '''
        Function simulate:
        Traces the rays through every component in the componentlist
        
        Inputs:
        None
        
        Outputs:
        The Rays that have been traced through all the components
        
        Notes:
        - The original Rays are modified by a call to this function
        - Calling this function will refresh the efficiency list and store the
            effiencies for each component
        '''
        # Get the source rays
        rays = self.source
        # Refresh the efficiency list
        self.effs = []
        
        # Iterate through each component
        for c in self.componentlist:
            
            # Trace the Rays through the component
            e = c.trace(rays,considerweights=self.weighting)
            
            # Some components return a list of efficiencies, check if this is the case
            if type(e) == list:
                self.effs.extend(e)
            else:
                self.effs.append(e)
            
            # Stop if there are no rays to continue going
            if (len(rays) == 0):
                break
        
        # Return the final Rays
        return rays
    
    
    ## Display Functions
    
    
    def displayEfficiency(self):
        '''
        Function displayEfficiency:
        Prints the efficiencies of the various components in an easily readable
            format
        
        Inputs:
        None
        Outputs:
        None
        
        Notes:
        - If there is no data in self.effs, this funcion will return without
            printing anything
        - The output consists of three parts:
            Method - How the photons were removed
            Local Eff - What percent of the photons that reached the component
                were removed
            Global Eff - What percent of the total photons were removed by this
                component specifically
        '''
        if (len(self.effs) == 0):
            return
        # Store the original number of rays
        l = self.effs[0][1]
        
        # Print the header
        print("Method:                        Local Percent:  Global Percent:")
        
        for e in self.effs:
            method = e[0].ljust(30)
            localeff = (e[1]-e[2])*100/e[1]
            globaleff = (e[1]-e[2])*100/l
            print("{0} {1:08.5f}%       {2:08.5f}%".format(method,localeff,globaleff))
        
        finalnum = self.effs[-1][2]
        
        print()
        print("Total Throughput: {0:08.5f}%".format(finalnum * 100 / l))
    
    
    
    
    
    
    
    
    