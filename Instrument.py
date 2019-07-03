import numpy as np
from prtp.Sources import Source
from prtp.Modification import Modification
from prtp.Combination import Combination
import astropy.units as u
import matplotlib.pyplot as plt

class Instrument:
    # The main class for simulating optical systems, each Instrument object
    # will contain a list of components you want to trace the rays to in order
    
    def __init__(self,source = None,weighting = False):
        '''
        Initialization Function:
        Creates the Instrument Objest
        
        Inputs:
        source - A Source Object, this will generate the Rays that you want to 
            send through the instrument. This can also be in Instrument object,
            if this is the case, this new Instrument will be given a source 
            that contains the Rays that the argument Instrument has generated.
        weighting - If False, the instrument will not weight photons
        
        Notes: 
        - componentlist is the list of all the components you will trace 
            photons to, in order.
        - effs is the list of efficiencies for each component, it will be filled
            when photons are traced to each component
        - self.rays will contain the result after calling simulate()
        '''
        if (type(source) == Instrument):
            if source.rays is None:
                raise Exception('Source Instrument has not yet been simulated')
            else:
                self.source = Source.sourceFromRays(source.rays)
        else:
            self.source = source
        self.weighting = weighting
        self.componentlist = []
        self.effs = []
        self.rays = None
    
    
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
    
    def addFocus(self,index=None):
        def focusfunc(rays,cw):
            rays.focusX()
        
        m = Modification(focusfunc)
        
        self.addComponent(m,index)
    
    def spectralResolution(self):
        '''
        Function spectralResolution:
        Returns the spectral resolution of the rays with dispersion in the
            x-direction, using the formula res=mean(x)/fwhm(x)
        
        Inputs:
        None
        
        Outputs:
        The spectral resolution as a float
        
        Notes:
        - There must be rays in the Instrument, meaning that it must have been
            simulated prior to the call to spectralResolution()
        '''
        if self.rays is None:
            raise ValueError('Instrument has not yet been simulated')
        else:
            return self.rays.spectralResolution()
    
    ## Misalignment Tests
    
    @u.quantity_input(min=u.mm,max=u.mm)
    def singleTranslateTest(self,index=0,min=-1*u.mm,max=1*u.mm,num=10,dim=1,plot=True,param='x'):
        '''
        Function singleTranslateTest:
        Performs repeated translational misalignments on the specified component
            and compares the FWHM after each one
        
        Inputs:
        index - The index in self.componentlist that refers to the component
            which you want to misalign
        min,max - This function will generate misalignment values equally spaced
            on the range [min,max]
        mun - The number of misalignment values you want to test
        dim - Can be 1, 2, or 3. Specifies if you want to translate in the 
            x-direction, y-direction, or z-direction, respectively.
        plot - boolean. If True, the misalignment and FWHM values will
            automatically be plotted in a scatter plot. If False, the values
            will just be returned
        param - The parameter in which you would like to plot the FWHM, can be
            x, y, z, l, m, n, ux, uy, or uz
        
        Outputs:
        misValues - The misalignment values that were tested, not as astropy
            Quantity but it gives values in mm
        resValues - The FWHM values produced by each misalignment value
        '''
        min = min.to(u.mm)
        max = max.to(u.mm)
        
        # Create an Instrument that simulates this up until the point of the
        # misalignment, this means all previous components will not have to be
        # simulated more than once, speeding up the test
        # prevInst = Instrument(self.source)
        # prevInst.componentlist = self.componentlist[:index]
        # prevInst.simulate()
        
        # misalignedInst = Instrument(prevInst)
        # misalignedInst.componentlist = self.componentlist[index:]        
        misalignedInst = Instrument(self.source)
        misalignedInst.componentlist = self.componentlist
        
        misValues = np.linspace(min.value,max.value,num)
        resValues = np.zeros(num)
        
        if type(misalignedInst.componentlist[index] == Combination):
        
            for i in range(num):
                # In each case, we must misalign the relevant component, simulate 
                # the instrument, then re-align the component
                if dim == 1:
                    Combination.translate(misalignedInst.componentlist[index],dx=misValues[i]*u.mm)
                    misalignedInst.simulate()
                    Combination.translate(misalignedInst.componentlist[index],dx=-1*misValues[i]*u.mm)
                elif dim == 2:
                    Combination.translate(misalignedInst.componentlist[index],dy=misValues[i]*u.mm)
                    misalignedInst.simulate()
                    Combination.translate(misalignedInst.componentlist[index],dy=-1*misValues[i]*u.mm)
                elif dim == 3:
                    Combination.translate(misalignedInst.componentlist[index],dz=-1*misValues[i]*u.mm)
                    misalignedInst.simulate()
                    Combination.translate(misalignedInst.componentlist[index],dz=misValues[i]*u.mm)
                else:
                    raise ValueError('dim must be 1, 2, or 3')
                resValues[i] = misalignedInst.getRays().fwhm(param=param)
            
        else:
            
            for i in range(num):
                # In each case, we must misalign the relevant component, simulate 
                # the instrument, then re-align the component
                if dim == 1:
                    misalignedInst.componentlist[index].translate(dx=misValues[i]*u.mm)
                    misalignedInst.simulate()
                    misalignedInst.componentlist[index].translate(dx=-1*misValues[i]*u.mm)
                elif dim == 2:
                    misalignedInst.componentlist[index].translate(dy=misValues[i]*u.mm)
                    misalignedInst.simulate()
                    misalignedInst.componentlist[index].translate(dy=-1*misValues[i]*u.mm)
                elif dim == 3:
                    misalignedInst.componentlist[index].translate(dz=misValues[i]*u.mm)
                    misalignedInst.simulate()
                    misalignedInst.componentlist[index].translate(dz=-1*misValues[i]*u.mm)
                else:
                    raise ValueError('dim must be 1, 2, or 3')
                resValues[i] = misalignedInst.getRays().fwhm(param=param)
        
        # Plot the Resolutions
        if plot:
            plt.figure()
            plt.scatter(misValues,resValues)
            plt.xlabel('Misalignment (mm)')
            plt.ylabel(param + ' FWHM')
            plt.show()
        
        return misValues,resValues
    
    
    @u.quantity_input(min=u.rad,max=u.rad)
    def singleUnitRotateTest(self,index=0,min=-1*u.deg,max=1*u.deg,num=10,axis=1,plot=True,param='x'):
        '''
        Function singleUnitRotateTest:
        Performs repeated unit-rotational misalignments on the specified component
            and compares the FWHM after each one
        
        Inputs:
        index - The index in self.componentlist that refers to the component
            which you want to misalign
        min,max - This function will generate misalignment values equally spaced
            on the range [min,max]
        mun - The number of misalignment values you want to test
        axis - Can be 1, 2, or 3. Specifies if you want to rotate about the 
            x-axis, y-axis, or z-axis, respectively.
        plot - boolean. If True, the misalignment and FWHM values will
            automatically be plotted in a scatter plot. If False, the values
            will just be returned
        param - The parameter in which you would like to plot the FWHM, can be
            x, y, z, l, m, n, ux, uy, or uz
        
        Outputs:
        misValues - The misalignment values that were tested, not as astropy
            Quantity but it gives values in degrees
        resValues - The FWHM values produced by each misalignment value
        '''
        min = min.to(u.rad)
        max = max.to(u.rad)
        
        # Create an Instrument that simulates this up until the point of the
        # misalignment, this means all previous components will not have to be
        # simulated more than once, speeding up the test
        # prevInst = Instrument(self.source)
        # prevInst.componentlist = self.componentlist[:index]
        # prevInst.simulate()
        
        # misalignedInst = Instrument(prevInst)
        # misalignedInst.componentlist = self.componentlist[index:]        
        misalignedInst = Instrument(self.source)
        misalignedInst.componentlist = self.componentlist
        
        misValues = np.linspace(min.value,max.value,num)
        resValues = np.zeros(num)
        
        if type(misalignedInst.componentlist[index] == Combination):
        
            for i in range(num):
                # In each case, we must misalign the relevant component, simulate 
                # the instrument, then re-align the component
                Combination.unitrotate(misalignedInst.componentlist[index],theta=misValues[i]*u.rad,axis=axis)
                misalignedInst.simulate()
                Combination.unitrotate(misalignedInst.componentlist[index],theta=-1*misValues[i]*u.rad,axis=axis)
                resValues[i] = misalignedInst.getRays().fwhm(param=param)
            
        else:
            
            for i in range(num):
                # In each case, we must misalign the relevant component, simulate 
                # the instrument, then re-align the component
                misalignedInst.componentlist[index].unitrotate(theta=misValues[i]*u.rad,axis=axis)
                misalignedInst.simulate()
                misalignedInst.componentlist[index].unitrotate(theta=-1*misValues[i]*u.rad,axis=axis)
                resValues[i] = misalignedInst.getRays().fwhm(param=param)
        
        # Plot the Resolutions
        if plot:
            plt.figure()
            plt.scatter((misValues*u.rad).to(u.deg),resValues)
            plt.xlabel('Misalignment (degrees)')
            plt.ylabel(param + ' FWHM')
            plt.show()
        
        return misValues,resValues
    
    
    @u.quantity_input(min=u.rad,max=u.rad)
    def singleRotateTest(self,index=0,min=-1*u.deg,max=1*u.deg,num=10,ux=1,uy=0,uz=0,plot=True,param='x'):
        '''
        Function singleRotateTest:
        Performs repeated unit-rotational misalignments on the specified component
            and compares the FWHM after each one
        
        Inputs:
        index - The index in self.componentlist that refers to the component
            which you want to misalign
        min,max - This function will generate misalignment values equally spaced
            on the range [min,max]
        mun - The number of misalignment values you want to test
        ux, uy, uz - These arguments describe the axis about which you want
            to rotate
        plot - boolean. If True, the misalignment and FWHM values will
            automatically be plotted in a scatter plot. If False, the values
            will just be returned
        param - The parameter in which you would like to plot the FWHM, can be
            x, y, z, l, m, n, ux, uy, or uz
        
        Outputs:
        misValues - The misalignment values that were tested, not as astropy
            Quantity but it gives values in degrees
        resValues - The FWHM values produced by each misalignment value
        '''
        min = min.to(u.rad)
        max = max.to(u.rad)
        
        # Create an Instrument that simulates this up until the point of the
        # misalignment, this means all previous components will not have to be
        # simulated more than once, speeding up the test
        # prevInst = Instrument(self.source)
        # prevInst.componentlist = self.componentlist[:index]
        # prevInst.simulate()
        
        # misalignedInst = Instrument(prevInst)
        # misalignedInst.componentlist = self.componentlist[index:]        
        misalignedInst = Instrument(self.source)
        misalignedInst.componentlist = self.componentlist
        
        misValues = np.linspace(min.value,max.value,num)
        resValues = np.zeros(num)
        
        if type(misalignedInst.componentlist[index] == Combination):
        
            for i in range(num):
                # In each case, we must misalign the relevant component, simulate 
                # the instrument, then re-align the component
                Combination.rotate(misalignedInst.componentlist[index],theta=misValues[i]*u.rad,ux=ux,uy=uy,uz=uz)
                misalignedInst.simulate()
                Combination.rotate(misalignedInst.componentlist[index],theta=-1*misValues[i]*u.rad,ux=ux,uy=uy,uz=uz)
                resValues[i] = misalignedInst.getRays().fwhm(param=param)
            
        else:
            
            for i in range(num):
                # In each case, we must misalign the relevant component, simulate 
                # the instrument, then re-align the component
                misalignedInst.componentlist[index].rotate(theta=misValues[i]*u.rad,ux=ux,uy=uy,uz=uz)
                misalignedInst.simulate()
                misalignedInst.componentlist[index].rotate(theta=-1*misValues[i]*u.rad,ux=ux,uy=uy,uz=uz)
                resValues[i] = misalignedInst.getRays().fwhm(param=param)
        
        # Plot the Resolutions
        if plot:
            plt.figure()
            plt.scatter((misValues*u.rad).to(u.deg),resValues)
            plt.xlabel('Misalignment (degrees)')
            plt.ylabel(param + ' FWHM')
            plt.show()
        
        return misValues,resValues
        
    
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
        rays = self.source.generateRays()
        
        if self.weighting:
            rays.addWeighting()

        # Refresh the efficiency list
        self.effs = []
        # Iterate through each component
        for c in self.componentlist:
            
            # Trace the Rays through the component
            e = c.trace(rays,considerweights=self.weighting)
            
            # Allow for some tracing to return None if Rays are never removed
            # by the component
            if e is None:
                continue
            
            # Some components return a list of efficiencies, check if this is the case
            if type(e) == list:
                self.effs.extend(e)
            else:
                self.effs.append(e)
            
            # Stop if there are no rays to continue going
            if (len(rays) == 0):
                break
        
        # Store and return the final Rays
        self.rays = rays
        return rays
    
    
    def getRays(self):
        '''
        Function getRays:
        Returns the rays after calling simulate.
        
        Inputs:
        None
        
        Outputs:
        Simulated Rays
        
        Notes:
        - Instrument.simulate() returns the rays, but self.rays also stores the
            rays in case the user wants to access them later
        '''
        return self.rays
    
    
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
        print("Method:                        Local Percent:   Global Percent:")
        
        for e in self.effs:
            method = e[0].ljust(30)
            if (e[1] == 0):
                localeff = 0.
                globaleff = 0.
            else:
                localeff = (e[1]-e[2])*100/e[1]
                globaleff = (e[1]-e[2])*100/l
            print("{0} {1:09.5f}%       {2:09.5f}%".format(method,localeff,globaleff))
        
        finalnum = self.effs[-1][2]
        
        print()
        print("Total Throughput: {0:08.5f}%".format(finalnum * 100 / l))
