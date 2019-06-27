import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import prtp.sources_helper as sources
import prtp.analyses as analyses
import prtp.surfacesf as surfacesf
import prtp.transformationsf as transformationsf
import prtp.woltsurf as wolt
from prtp.analyses import analyticYPlane,analyticXPlane,analyticImagePlane

class Rays:
    
    
    ## Creation Functions:
    # These functions create Rays Objects
    
    
    def __init__(self,x=[],y=[],z=[],l=[],m=[],n=[],ux=[],uy=[],uz=[],
                weighting=False,wave=None,order=None):
        '''
        Function __init__:
        Initializes a Rays Object
        
        Inputs:
        x,y,z - The spatial location of the photons
        l,m,n - The direction components of the photons
        ux,uy,uz - The components of the normal vector of the most recent 
            surface these photons were traced to
        weighting - Used by Instrument Objects, see notes
        wave - The wavelength of the photons, can be a singular value or a
            numpy array that has the same length as the Rays object
        order - The orders of the photons, can be a singular value or a
            numpy array that has the same length as the Rays object
        
        Notes:
        - If weighting is True, each photon will immediately be given a weight
            of 1. Objects like Gratings and Detectors have reflectivities and
            quantum efficiencies (respectively) that will remove photons based
            on a probability. With weighted photons, no photons will be removed
            by these objects but their weights will be adjusted. A photon's 
            weight is simply the chance that this specific photon made it this
            far through an instrument. This allows us to get accurate instrument
            efficiency values using relatively fewer photons.
        - If wave or order are singular values, that value will be given to
            every photon in the array.
        '''
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.l = np.array(l)
        self.m = np.array(m)
        self.n = np.array(n)
        self.ux = np.array(ux)
        self.uy = np.array(uy)
        self.uz = np.array(uz)
        
        self.weighting = weighting
        if weighting:
            self.weights = np.ones(len(self.x))
        else:
            self.weights = None
        
        # This call to the __len__ method makes sure the inputs were all of the
        # same length
        self.__len__()
        self.tags = []
        self.params = []
        self.startingeff = None
        
        
        if (wave is not None):
            if type(wave) == np.ndarray:
                self.wave = wave
            else:
                self.wave = np.ones(len(self.x)) * wave
        else:
            self.wave = None
        if (order is not None):
            if type(self.order) == np.ndarray:
                self.order = order
            else:
                self.order = np.ones(len(self.x)) * order
        else:
            self.order = None
    
    
    @classmethod
    def pointsource(cls,ang,num):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.pointsource(
                                           ang,num)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def circularbeam(cls,rad,num):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.circularbeam(
                                           rad,num)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def annulus(cls,rin,rout,num,zhat=-1.):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.annulus(
                                           rin,rout,num,zhat)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def subannulus(cls,rin,rout,dphi,num,zhat=1.):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.subannulus(
                                           rin,rout,dphi,num,zhat)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def xslit(cls,xin,xout,num,zhat=-1.):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.xslit(
                                           xin,xout,num,zhat)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def rectarray(cls,xsize,ysize,num):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.rectArray(
                                           xsize,ysize,num)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def convergingbeam(cls,zset,rin,rout,tmin,tmax,num,lscat):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.convergingbeam(
                                           zset,rin,rout,tmin,tmax,num,lscat)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def convergingbeam2(cls,zset,xmin,xmax,ymin,ymax,num,lscat):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.convergingbeam2(
                                           zset,xmin,xmax,ymin,ymax,num,lscat)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def rectbeam(cls,xhalfwidth,yhalfwidth,num):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.rectbeam(
                                           xhalfwidth,yhalfwidth,num)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    def copy(self):
        '''
        Function copy:
        Returns a copy of the Rays object
        
        Inputs:
        Nothing
        
        Outputs:
        An identical Rays object
        
        Notes:
        -The returned Rays object is a deep copy of the original, modifying one
        should have no effect on the other
        '''
        x = np.copy(self.x)
        y = np.copy(self.y)
        z = np.copy(self.z)
        l = np.copy(self.l)
        m = np.copy(self.m)
        n = np.copy(self.n)
        ux = np.copy(self.ux)
        uy = np.copy(self.uy)
        uz = np.copy(self.uz)
        new_rays = Rays(x,y,z,l,m,n,ux,uy,uz)
        new_rays.tags = self.tags.copy()
        new_rays.params = self.params.copy()
        if self.wave is not None:
            new_rays.wave = self.wave.copy()
        if self.order is not None:
            new_rays.order = self.order.copy()
        return new_rays
    
    def makecopy(self,rays):
        '''
        Function makecopy:
        Makes this Rays object a copy of another Rays object. This is a very
        niche function that is needed when tracing through Combination
        objects.
        
        Inputs:
        rays - The rays object that you want to make self a copy of
        
        Outputs:
        None
        '''
        self.x = np.copy(rays.x)
        self.y = np.copy(rays.y)
        self.z = np.copy(rays.z)
        self.l = np.copy(rays.l)
        self.m = np.copy(rays.m)
        self.n = np.copy(rays.n)
        self.ux = np.copy(rays.ux)
        self.uy = np.copy(rays.uy)
        self.uz = np.copy(rays.uz)
        self.tags = rays.tags.copy()
        self.params = rays.params.copy()
        self.weighting = rays.weighting
        if self.weighting:
            self.weights = rays.weights
        else:
            self.weights = None
        self.wave = rays.wave
        self.order = rays.order
    
    def pullrays(self,rays,tarray):
        '''
        Function pullrays:
        Takes photons from another Rays object and copies them into this one
        
        Inputs:
        rays - The rays object from which you want to pull photons
        tarray - A trutharray that dictates which photons are to be pulled.
            For example, if trutharray[0] is True, the first photon in self
            will be changed into a copy of the first photon in rays. If
            trutharray[1] is False, the second photon in self will be kept
            as it was.
        
        Outputs:
        None
        
        Notes:
        - self and rays must contain the same number of photons.
        - The function was created for use in GratingStack's trace functions
        - This might not work well for Rays objects with Tags or Params
        '''
        self.x[tarray] = np.copy(rays.x[tarray])
        self.y[tarray] = np.copy(rays.y[tarray])
        self.z[tarray] = np.copy(rays.z[tarray])
        self.l[tarray] = np.copy(rays.l[tarray])
        self.m[tarray] = np.copy(rays.m[tarray])
        self.n[tarray] = np.copy(rays.n[tarray])
        self.ux[tarray] = np.copy(rays.ux[tarray])
        self.uy[tarray] = np.copy(rays.uy[tarray])
        self.uz[tarray] = np.copy(rays.uz[tarray])
        
        for i in range(len(self.tags)):
            value = self.tags[i][1]
            value[tarray] = rays.getTag(self.tags[i][0])[tarray]
            self.tags[i] = (self.tags[i][0],value)
        for i in range(len(self.params)):
            value = self.params[i][1]
            value[tarray] = rays.getTag(self.params[i][0])[tarray]
            self.params[i] = (self.params[i][0],value)
        
        if self.weighting:
            self.weights[tarray] = rays.weights[tarray]
        
        if self.wave is not None:
            self.wave[tarray] = rays.wave[tarray]
        if self.order is not None:
            self.order[tarray] = rays.order[tarray]
    
    def split(self,trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function copy:
        Returns a new Rays object containing the photons which satisy some
        condition
        
        Inputs:
        trutharray,tags,delim,orcombination - Arguments needed to develop a trutharray which specifies the photons which should be copied.
        See documentation for combineTags to see how it works
        
        Outputs:
        A Rays object containing the photons which satisy trutharray
        
        Notes:
        -The returned Rays object is a deep copy of the original, modifying one
        should have no effect on the other
        -The original Rays object is unaffected by a call to split()
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        
        # Copy main attributes
        x = np.copy(self.x[tarray])
        y = np.copy(self.y[tarray])
        z = np.copy(self.z[tarray])
        l = np.copy(self.l[tarray])
        m = np.copy(self.m[tarray])
        n = np.copy(self.n[tarray])
        ux = np.copy(self.ux[tarray])
        uy = np.copy(self.uy[tarray])
        uz = np.copy(self.uz[tarray])
        new_rays = Rays(x,y,z,l,m,n,ux,uy,uz)
        
        # Copy Tags and Parameters
        for tag in self.tags:
            new_rays.addTag(tag[0],tag[1][tarray])
        for param in self.params:
            new_rays.addParam(param[0],param[1][tarray])
        
        # Consider Weighting
        if self.weighting:
            new_rays.weighting = True
            new_rays.weights = np.copy(self.weights[tarray])
        
        if self.wave is not None:
            new_rays.wave = np.copy(self.wave[tarray])
        if self.order is not None:
            new_rays.order = np.copy(self.order[tarray])
        
        return new_rays
    
    ## Info Functions:
    # These functions give information about the Rays Object
    
    def __len__(self):
        '''
        Function __len__:
        Returns the length of the x,y,z,l,m,n,ux,uy,uz parameters for this 
        Rays object
        
        Inputs:
        Nothing
        
        Outputs:
        length - The length of the parameter arrays, if they are of the same 
        length
        
        Notes:
        -If the x,y,z,l,m,n,ux,uy,uz arrays are not all of the same length,
        this function raises an exception
        -This function is called using the syntax:
        >>len(r)
        where r is a Rays object
        '''
        if (len(self.x) == len(self.y) == len(self.z) == len(self.l) == 
            len(self.m) == len(self.n) == len(self.ux) == len(self.uy) ==
            len(self.uz)):
            return len(self.x)
        else:
            raise Exception('The x,y,x,l,m,n,ux,uy,uz Parameters must be of the same length!')
    
    def length(self,considerweights=False):
        '''
        Function length:
        Gets the length of this Rays Object. Differs from __len__ in that this
        can consider the weighting of the Rays Object
        
        Inputs:
        considerweights - If True, this function will return the sum of the 
            photons' weights. Otherwise, behaves identically to __len__
        
        Outputs:
        None
        '''
        if (considerweights):
            if self.weights is None:
                raise Exception('This Rays Object has no weighting')
            else:
                return np.sum(self.weights)
        else:
            return len(self)
    
    
    ## Weighting Functions
    # Functions that deal with photons' weights
    
    def addWeighting(self):
        '''
        Function addWeighting:
        Adds weights to a Rays object even if it was not initialized with weights
        
        Inputs:
        None
        
        Outputs:
        None
        
        Notes:
        If this function is called on a Rays Object that already has weighting,
            the existing weights will be reset (back to 1. for each photon).
        '''
        self.weighting = True
        self.weights = np.ones(len(self.x))
    
    
    def removeWeighting(self):
        '''
        Function removeWeighting:
        Removes weights from a Rays object
        
        Inputs:
        None
        
        Outputs:
        None
        
        Notes:
        This function will not do anything if the Rays object did not have 
            weights to begin with.
        '''
        self.weighting = False
        self.weights = None
    
    
    ## Analysis Functions:
    # These functions are calls to prtp.analyses
    
    def centroid(self,weights=None):
        return analyses.centroid(self,weights)
    
    
    def rmscentroid(self,weights=None):
        return analyses.rmscentroid(self,weights)
    
    
    def rmsX(self,weights=None):
        return analyses.rmscentroid(self,weights)
    
    
    def rmsY(self,weights=None):
        return analyses.rmsX(self,weights)
    
    
    def rho(self,weights=None,cent=False):
        return analyses.rho(self,weights,cent)
    
    
    def rhocdf(self,weights=None,cent=False):
        return analyses.rhocdf(self,weights,cent)
    
    
    def hpd(self,weights=None):
        return analyses.hpd(self,weights)
    
    
    def hpdY(self,weights=None):
        return analyses.hpdY(self,weights)
    
    
    def analyticImagePlane(self,weights=None):
        return analyses.analyticImagePlane(self,weights)
    
    
    def analyticXPlane(self,weights=None):
        return analyses.analyticXPlane(self,weights)
    
    
    def analyticYPlane(self,weights=None):
        return analyses.analyticYPlane(self,weights)
    
    
    def indAngle(self,ind=None,normal=None):
        return analyses.indAngle(self,ind,normal)
    
    
    def grazeAngle(self,ind=None):
        return analyses.grazeAngle(self,ind)
    
    
    def interpolateVec(self,I,Nx,Ny,xr=None,yr=None,method='linear',polar=False):
        return analyses.interpolateVec(self,I,Nx,Ny,xr,yr,method,polar)
    
    # wavefront requires the no longer used utilities.imaging.man package
    # def wavefront(self,Nx,Ny,method='cubic',polar=False):
    #     return analyses.wavefront(self,Nx,Ny,method,polar)
    
    
    # measurePower requires the no longer used utilities.imaging.man package
    # def measurePower(self,Nx,Ny,method='linear'):
    #     return analyses.measurePower(self,Nx,Ny,method)
    
    
    ## Surfaces Functions:
    # These functions call to prtp.surfacesf
    
    def flat(self):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.flat(self)
    
    
    def sphere(self,rad,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.sphere(self,rad)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def cyl(self,rad,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.cyl(self,rad)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))

    def cylconic(self,rad,k,eliminate='nan',maxiter=12):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.cylconic(self,rad,k,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def conic(self,rad,k,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.conic(self,rad,k)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def paraxial(self,f):
        if (f==0):
            raise ValueError('f cannot be 0')
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.paraxial(self,f)
    
    def paraxialY(self,f):
        if (f==0):
            raise ValueError('f cannot be 0')
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.paraxial(self,f)
    
    def torus(self,rin,rout,eliminate='nan',maxiter=12):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.torus(self,rin,rout,eliminate,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def conicplus(self,r,k,n,p,Np,eliminate='nan',maxiter=12):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = surfacesf.conicplus(self,r,k,n,p,Np,eliminate,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def focus(self,fn,weights=None):
        return surfacesf.focus(self,fn,weights)
    
    def focusX(self,weights=None):
        return self.focus(analyticXPlane,weights)
    
    def focusY(self,weights=None):
        return self.focus(analyticYPlane,weights)
    
    def focusI(self,weights=None):
        return self.focus(analyticImagePlane,weights)
    

    ## Basic Motion Functions:
    
    def translate(self,dx=0,dy=0,dz=0, trutharray=None,tags=None,delim=None,orcombination=True,coords=None):
        '''
        Function translate:
        Moves the position of a Rays object by a specified amount
        
        Inputs:
        dx - The amount to change in x
        dy - The amount to change in y
        dz - The amount to change in z
        
        Outputs:
        Nothing
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        self.x[tarray] += dx
        self.y[tarray] += dy
        self.z[tarray] += dz
        
        if coords is not None:
            #Define rotation and translation matrices
            tranm = transformationsf.translationM(-dx,-dy,-dz)
            tranmi = transformationsf.translationM(dx,dy,dz)
            #Dot rotation into forward transform
            coords[1] = np.dot(tranm,coords[1])
            coords[3] = np.dot(coords[3],tranmi)
    
    
    def rotate(self,dl=0,dm=0,dn=0, trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function rotate:
        Moves the direction of a Rays object by a specified amount
        
        Inputs:
        dl - The amount to change in l
        dm - The amount to change in m
        dn - The amount to change in n
        
        Outputs:
        Nothing
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        self.l[tarray] += dl
        self.m[tarray] += dm
        self.n[tarray] += dn
    
    
    def rotatenormal(self,dux=0,duy=0,duz=0, trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function rotatenormal:
        Moves the normal vector of a Rays object by a specified amount
        
        Inputs:
        dux - The amount to change in ux
        duy - The amount to change in uy
        duz - The amount to change in uz
        
        Outputs:
        Nothing
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        self.ux += dux
        self.uy += duy
        self.uz += duz
    
    def set(self,x=None,y=None,z=None,
                  l=None,m=None,n=None,
                  ux=None,uy=None,uz=None, trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function set:
        Changes the components of a Rays object to be a new array
        This function "sets the position of each photon"
        
        Inputs:
        x,y,z,l,m,n,ux,uy,uz - New arrays to take the place of the current x,y,z,l,m,n,ux,uy, or uz
        
        Outputs:
        Nothing
        
        Notes:
        - This function differs from move. This function takes in new arrays while Move takes in constants and makes uniform arrays
        - This function allows the user to set specific values for each photon without accessing the class variables
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        size = len(tarray)
        if x is not None:
            if (len(x) == size):
                self.x[tarray] = x
            else:
                raise ValueError('x input not of the correct size')
        if y is not None:
            if (len(y) == size):
                self.y[tarray] = y
            else:
                raise ValueError('y input not of the correct size')
        if z is not None:
            if (len(z) == size):
                self.z[tarray] = z
            else:
                raise ValueError('z input not of the correct size')
        if l is not None:
            if (len(l) == size):
                self.l[tarray] = l
            else:
                raise ValueError('l input not of the correct size')
        if m is not None:
            if (len(m) == size):
                self.m[tarray] = m
            else:
                raise ValueError('m input not of the correct size')
        if n is not None:
            if (len(n) == size):
                self.n[tarray] = n
            else:
                raise ValueError('n input not of the correct size')
        if ux is not None:
            if (len(ux) == size):
                self.ux[tarray] = ux
            else:
                raise ValueError('ux input not of the correct size')
        if uy is not None:
            if (len(uy) == size):
                self.uy[tarray] = uy
            else:
                raise ValueError('uy input not of the correct size')
        if uz is not None:
            if (len(uz) == size):
                self.uz[tarray] = uz
            else:
                raise ValueError('uz input not of the correct size')
    
    
    def move(self,x=None,y=None,z=None,
                  l=None,m=None,n=None,
                  ux=None,uy=None,uz=None, trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function move:
        Changes the components of a Rays object to be a new array
        This function "Moves all Rays to a new position"
        
        Inputs:
        x,y,z,l,m,n,ux,uy,uz - Constants that will specify the new values of each parameter
        
        Outputs:
        Nothing
        
        Notes:
        - This function differs from set. This function takes in constants while Set takes in new arrays
        - For example, if the x argument is set to 3., every photon in the Rays object will be given an x-position of 3.
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        size = len(tarray)
        if x is not None:
            self.x[tarray] = np.ones(size) * x
        if y is not None:
            self.y[tarray] = np.ones(size) * y
        if z is not None:
            self.z[tarray] = np.ones(size) * z
        if l is not None:
            self.l[tarray] = np.ones(size) * l
        if m is not None:
            self.m[tarray] = np.ones(size) * m
        if n is not None:
            self.n[tarray] = np.ones(size) * n
        if ux is not None:
            self.ux[tarray] = np.ones(size) * ux
        if uy is not None:
            self.uy[tarray] = np.ones(size) * uy
        if uz is not None:
            self.uz[tarray] = np.ones(size) * uz
    

    ## Transformation Functions:
    # These functions call to prtp.transformationsf
    
    def rotatevector(self,type,theta,axis):
        '''
        Function rotatevector:
        Rotates a part of the Rays object about a unit-axis
        
        Inputs:
        type - A string describing which part of the Rays object you want to rotate. Options are listed below (case-insensitive):
            pos - The x,y,z positions are rotated
            posn - The x,y,z positions are rotated
            position - The x,y,z positions are rotated
            dir - The l,m,n directions are rotated
            direction - The l,m,n directions are rotated
            norm - The ux,uy,uz directions are rotated
            normal - The ux,uy,uz directions are rotated
            all - x,y,z,l,m,n,ux,uy,uz are all rotated
        theta - The angle through which you want to rotate (in radians)
        axis - The integer 1,2 or 3 corresponding to the axis about which you want to rotate. 1,2,3 correspond to x,y,z, respectively
        
        Outputs:
        None
        
        Notes:
        - theta can be an array with the same length as the Rays object. If it is, each photon will be rotated a unique amount defined by the theta array
        '''
        type = type.lower()
        if (type == 'pos' or type == 'posn' or type == 'position'):
            self.x,self.y,self.z = transformationsf.rotatevector(self.x,self.y,self.z,theta,axis)
        
        elif (type == 'dir' or type == 'direction'):
            self.l,self.m,self.n = transformationsf.rotatevector(self.l,self.m,self.n,theta,axis)
        
        elif (type == 'norm' or type == 'normal'):
            self.ux,self.uy,self.uz = transformationsf.rotatevector(self.ux,self.uy,self.uz,theta,axis)
        
        elif (type == 'all'):
            self.x,self.y,self.z = transformationsf.rotatevector(self.x,self.y,self.z,theta,axis)
            self.l,self.m,self.n = transformationsf.rotatevector(self.l,self.m,self.n,theta,axis)
            self.ux,self.uy,self.uz = transformationsf.rotatevector(self.ux,self.uy,self.uz,theta,axis)
        
        else:
            raise ValueError('type argument is not a valid value. Accepted values include pos,posn,position,dir,direction,norm,normal, and all')

    
    def rotatecustomaxis(self,type,theta,ax,ay,az):
        '''
        Function rotatecustomaxis:
        Rotates a part of the Rays object about a user-defined axis
        
        Inputs:
        type - A string describing which part of the Rays object you want to rotate. Options are listed below (case-insensitive):
            pos - The x,y,z positions are rotated
            posn - The x,y,z positions are rotated
            position - The x,y,z positions are rotated
            dir - The l,m,n directions are rotated
            direction - The l,m,n directions are rotated
            norm - The ux,uy,uz directions are rotated
            normal - The ux,uy,uz directions are rotated
            all - x,y,z,l,m,n,ux,uy,uz are all rotated
        theta - The angle through which you want to rotate (in radians)
        ax,ay,az - The x,y,z components of the axis about which you want to rotate
        
        Outputs:
        None
        
        Notes:
        - ax, ay, and az can be floats or arrays with the same length as the Rays object. If they are arrays, each photon will be rotated about a unique axis
        - theta can also be an array with the same length as the Rays object to rotate each photon a unique amount
        '''
        type = type.lower()
        if (type == 'pos' or type == 'posn' or type == 'position'):
            self.x,self.y,self.z,ax,ay,az = transformationsf.rotateaxis(self.x,self.y,self.z,ax,ay,az,theta)
        
        elif (type == 'dir' or type == 'direction'):
            self.l,self.m,self.n,ax,ay,az = transformationsf.rotateaxis(self.l,self.m,self.n,ax,ay,az,theta)
        
        elif (type == 'norm' or type == 'normal'):
            self.ux,self.uy,self.uz,ax,ay,az = transformationsf.rotateaxis(self.ux,self.uy,self.uz,ax,ay,az,theta)
        
        elif (type == 'all'):
            self.x,self.y,self.z,ax,ay,az = transformationsf.rotateaxis(self.x,self.y,self.z,ax,ay,az,theta)
            self.l,self.m,self.n,ax,ay,az = transformationsf.rotateaxis(self.l,self.m,self.n,ax,ay,az,theta)
            self.ux,self.uy,self.uz,ax,ay,az = transformationsf.rotateaxis(self.ux,self.uy,self.uz,ax,ay,az,theta)
        
        else:
            raise ValueError('type argument is not a valid value. Accepted values include pos,posn,position,dir,direction,norm,normal, and all')

    
    def reflect(self,trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function reflect:
        Reflects the rays off of their current surface
        
        Inputs:
        Used if only some photons are to be reflected, see 
        documentation for combinetags
        
        Outputs:
        Nothing
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        l,m,n = transformationsf.reflect(
            self.l[tarray],self.m[tarray],self.n[tarray], 
            self.ux[tarray],self.uy[tarray],self.uz[tarray])
        self.l[tarray] = l
        self.m[tarray] = m
        self.n[tarray] = n

    
    def refract(self,n1,n2, trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function reflect:
        Reflects the rays off of their current surface
        
        Inputs:
        n1 - The index refraction of the current medium
        n2 - The index of refraction of the medium the rays are going into
        Other inputs - Used if only some photons are to be reflected, see 
        documentation for combinetags
        
        Outputs:
        Nothing
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        l,m,n,ux,uy,uz = transformationsf.refract(     self.l[tarray],self.m[tarray],self.n[tarray], self.ux[tarray],self.uy[tarray],self.uz[tarray],n1,n2)
        self.l[tarray] = l
        self.m[tarray] = m
        self.n[tarray] = n
        self.ux[tarray] = ux
        self.uy[tarray] = uy
        self.uz[tarray] = uz

    
    def transform(self,tx,ty,tz,rx,ry,rz,coords=None, trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function transform:
        Performs a Coordinate Transformation on the Rays
        
        Inputs:
        tx,ty,tz - The translations of x, y, and z
        rx,ry,rz - The rotations around l, m, and n
        coords - If specified, updates a rotation matrix
        other inputs - Used if only some photons are to be transformed, see 
        documentation for combinetags
        
        Outputs:
        Nothing
        
        Notes:
        -Translations are performed first
        -If the user only desires to move the the rays' positions, 
        rays.translate is much faster
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        x,y,z,l,m,n,ux,uy,uz = transformationsf.transform(self.x[tarray],self.y[tarray],self.z[tarray],     self.l[tarray],self.m[tarray],self.n[tarray], self.ux[tarray],self.uy[tarray],self.uz[tarray], tx,ty,tz,rx,ry,rz,coords)
        self.x[tarray] = x
        self.y[tarray] = y
        self.z[tarray] = z
        self.l[tarray] = l
        self.m[tarray] = m
        self.n[tarray] = n
        self.ux[tarray] = ux
        self.uy[tarray] = uy
        self.uz[tarray] = uz


    def itransform(self,tx,ty,tz,rx,ry,rz,coords=None, trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function itransform:
        Performs an Inverse Coordinate Transformation on the Rays
        
        Inputs:
        tx,ty,tz - The translations of x, y, and z
        rx,ry,rz - The rotations around l, m, and n
        coords - If specified, updates a rotation matrix
        other inputs - Used if only some photons are to be transformed, see 
        documentation for combinetags
        
        Outputs:
        Nothing
        
        Notes:
        -Behaves as the opposite of transformation, that is:
        >>rays.transform(tx,ty,tz,rx,ry,rz)
        >>rays.itransform(tx,ty,tz,rx,ry,rz)
        will have no net effect on the rays (itransform undoes transform)
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        x,y,z,l,m,n,ux,uy,uz = transformationsf.itransform(self.x[tarray],self.y[tarray],self.z[tarray],     self.l[tarray],self.m[tarray],self.n[tarray], self.ux[tarray],self.uy[tarray],self.uz[tarray], tx,ty,tz,rx,ry,rz,coords)
        self.x[tarray] = x
        self.y[tarray] = y
        self.z[tarray] = z
        self.l[tarray] = l
        self.m[tarray] = m
        self.n[tarray] = n
        self.ux[tarray] = ux
        self.uy[tarray] = uy
        self.uz[tarray] = uz


    def radgrat(self,dpermm,order,wave,eliminate='nan', trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function radgrat:
        Sends the rays through a radial grating
        
        Inputs:
        dpermm - The groove density of the grating
        order - The order of the reflected photon
        wave - The wavelength of the incoming photons
        other inputs - Used if only some photons are to be transformed, see 
        documentation for combinetags
        
        Outputs:
        Nothing
        
        Notes:
        -Any of the inputs (dpermm,order, or wave) can be arrays, if an array
        is passed into the function, it must be the same length as the array
        of photons going onto the grating. Each element in the wave array gives
        the wavelength of an individual photon (same for dpermm and order).
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        
        # Check which optional arguments are numpy arrays:
        if (type(dpermm) == np.ndarray):
            dpermm = dpermm[tarray]
        if (type(order) == np.ndarray):
            order = order[tarray]
        if (type(wave) == np.ndarray):
            wave = wave[tarray]
            
        x,y,l,m,n = transformationsf.radgrat(self.x[tarray],self.y[tarray],self.l[tarray],self.m[tarray],self.n[tarray],dpermm,order,wave)
        self.x[tarray] = x
        self.y[tarray] = y
        self.l[tarray] = l
        self.m[tarray] = m
        self.n[tarray] = n
        
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))


    def grat(self,d,order,wave,eliminate='nan', trutharray=None,tags=None,delim=None,orcombination=True):
        '''
        Function grat:
        Sends the rays through a parallel grating
        
        Inputs:
        d - The groove density of the grating
        order - The order of the reflected photon
        wave - The wavelength of the incoming photons
        other inputs - Used if only some photons are to be transformed, see 
        documentation for combinetags
        
        Outputs:
        Nothing
        
        Notes:
        -Any of the inputs (d,order, or wave) can be arrays, if an array
        is passed into the function, it must be the same length as the array
        of photons going onto the grating. Each element in the wave array gives
        the wavelength of an individual photon (same for d and order).
        '''
        if tags is not None:
            # Combine Tags
            tarray = self.combineTags(tags,delim,orcombination)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            tarray = (trutharray != 0)
        else:
            tarray = np.full((len(self)), True, dtype=bool)
        
        # Check which optional arguments are numpy arrays:
        if (type(d) == np.ndarray):
            d = d[tarray]
        if (type(order) == np.ndarray):
            order = order[tarray]
        if (type(wave) == np.ndarray):
            wave = wave[tarray]
        
        l,m,n = transformationsf.grat(self.l[tarray],self.m[tarray],self.n[tarray],d,order,wave)
        self.l[tarray] = l
        self.m[tarray] = m
        self.n[tarray] = n
        
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))


    def applyT(self,coords,inverse=False):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = transformationsf.applyT(self,coords,inverse)
    
    
    def beckmann_scatter(self, h, rho, ripple, nind=1., lam=None):
        '''
        Add Beckmann scatter to a set of rays a la the SCATTER function
        in Webster Cash's IRT.
        '''
        if len(self) == 0:
            return
    
        # Define surface normal components.
        nx = self.ux
        ny = self.uy
        nz = self.uz
        n = np.array([nx, ny, nz]).transpose()
    
        # Define direction cosines.
        qx = self.l
        qy = self.m
        qz = self.n
        q = np.array([qx, qy, qz]).transpose()
    
        # Compute some cross products.
        op = np.cross(n, q)
        
        norm = np.linalg.norm(op,axis=1)
        op[:,0] /= norm
        op[:,1] /= norm
        op[:,2] /= norm
        
        ip = np.cross(q, op)
        
        norm = np.linalg.norm(ip,axis=1)
        ip[:,0] /= norm
        ip[:,1] /= norm
        ip[:,2] /= norm
    
        dot = np.sum((q*n),1)
        
        nel = len(self)
        
        dq = np.zeros(np.shape(q))
        
        # Create *h* array.
        if type(h) is int or float:
            h = [h]
            nh = 1
        else:
            nh = len(list(h))
    
        h = np.array(h)
    
        # Create *rho* array.
        if type(rho) is int or float:
            rho = [rho]
    
        rho = np.array(rho)
    
        for ih in range(nh):
            if h[ih] <= 0.:
                continue
            hran = np.random.rand(nel)
            lamn = lam / nind
            gsc = 4. * np.pi * h[ih] * dot / lamn
            gsc = np.exp(-gsc * gsc)
            scat = [hran < gsc]
    
            grad = ((1. / np.random.rand(nel)) - 1.) * lamn/rho[ih]
    
            gthet = np.random.rand(nel) * 2 * np.pi
            dthet = grad * np.sin(gthet) / dot
            dphi = grad * np.cos(gthet)
    
            dq = dq + (dthet*ip + dphi*op) * scat
    
        ripr = np.random.randn(nel) * ripple
        rth = np.random.rand(nel) * 2 * np.pi
        dth = ripr * np.sin(rth)
        dph = ripr * np.cos(rth) * dot
        
        dq[:,0] += (dth*ip[:,0]) + (dph*op[:,0])        
        dq[:,1] += (dth*ip[:,1]) + (dph*op[:,1])        
        dq[:,2] += (dth*ip[:,2]) + (dph*op[:,2])
        
        q = q + dq
        
        norm = np.linalg.norm(q,axis=1)
        q[:,0] /= norm
        q[:,1] /= norm
        q[:,2] /= norm
        
        q = q.transpose()
    
        self.l = np.array(q[0], order='F')
        self.m = np.array(q[1], order='F')
        self.n = np.array(q[2], order='F')
    
        return
    
    
    ## Wolter Surfaces Functions:
    # These functions call to prtp.woltsurf
    
    def wolterprimary(self,r0,z0,psi=1.,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.wolterprimary(self,r0,z0,psi,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def woltersecondary(self,r0,z0,psi=1.,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.woltersecondary(self,r0,z0,psi,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def woltersine(self,r0,z0,amp,freq,maxiter=12,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.woltersine(self,r0,z0,amp,freq,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def wolterprimLL(self,r0,z0,zmax,zmin,dphi,coeff,axial,az,cnum,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.wolterprimLL(self,r0,z0,zmax,zmin,dphi,coeff,axial,az,cnum,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def woltersecLL(self,r0,z0,psi,zmax,zmin,dphi,coeff,axial,az,cnum,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.woltersecLL(self,r0,z0,zmax,zmin,dphi,coeff,axial,az,cnum,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def wsprimary(self,alpha,z0,psi,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.wsprimary(self,alpha,z0,psi,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def wssecondary(self,alpha,z0,psi,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.wssecondary(self,alpha,z0,psi,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def spoCone(self,R0,tg,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.spoCone(self,R0,tg)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def ellipsoidWoltLL(self,r0,z0,psi,S,zmax,zmin,dphi,coeff,axial,ax,cnum,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.ellipsoidWoltLL(self,r0,z0,psi,S,zmax,zmin,dphi,coeff,axial,ax,cnum,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def wsprimaryBack(self,alpha,z0,psi,thick,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.wsprimaryBack(self,alpha,z0,psi,thick,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    def wssecondaryBack(self,alpha,z0,psi,thick,maxiter=10,eliminate='nan'):
        self.x,self.y,self.z,self.l,self.m,self.n,self.ux,self.uy,self.uz = wolt.wssecondaryBack(self,alpha,z0,psi,thick,maxiter)
        if eliminate.lower() == 'remove':
            self.remove(np.isnan(self.x))
    
    
    ## Update Functions:
    # These functions update the Rays object, this section does NOT include 
    # motion functions or functions that pass photons through optical systems.
    
    def remove(self,trutharray=None,tags=None,delim=None):
        '''
        Function remove:
        Removes photons from a Rays object given a condition
        
        Inputs:
        trutharray - The array that specifies which photons will be removed, a 1 (True) will be removed while a 0 (False) will be saved
        tags, delim - If these arguments are specified, a trutharray will be
        generated using the self.combinetags(tags,delim) method, this trutharray
        will be used in the same way the previous argument would be
        
        Outputs:
        Nothing
        
        Notes:
        -If the user tries to specify a trutharray, tags, and a delim, only the
        trutharray will be used. There are ways to combine a trutharray with
        existing tags (like creating a new tag), but that it not handled by this
        function
        '''
        if tags is not None:
            trutharray = self.combineTags(tags,delim)
        elif trutharray is not None:
            # Convert trutharray from a numerical to a boolean array
            trutharray = (trutharray != 0)
        else:
            return
        
        # Switch so that a true will be removed while a false will be saved
        trutharray = np.logical_not(trutharray)
            
        self.x = self.x[trutharray]
        self.y = self.y[trutharray]
        self.z = self.z[trutharray]
        self.l = self.l[trutharray]
        self.m = self.m[trutharray]
        self.n = self.n[trutharray]
        self.ux = self.ux[trutharray]
        self.uy = self.uy[trutharray]
        self.uz = self.uz[trutharray]
        # Iterate through tags and create new tags according to trutharray
        for i in range(len(self.tags)):
            tag = self.tags[i]
            tag = (tag[0],tag[1][trutharray])
            self.tags[i] = tag
        # Do the same for params
        for i in range(len(self.params)):
            param = self.params[i]
            param = (param[0],param[1][trutharray])
            self.params[i] = param
        # Remove from weighting, if necessary
        if (self.weighting):
            self.weights = self.weights[trutharray]
        
        if self.wave is not None:
            self.wave = self.wave[trutharray]
        if self.order is not None:
            self.order = self.order[trutharray]
    
    
    def probRemove(self,probability=1.,considerweights=False):
        '''
        Function probRemove:
        Removes photons from a Rays object based on a probability. It is the
        probability of survival, so a probability of 1 means the photon will 
        always remain and a probability of 0 means that the photon will always
        be removed.
        
        Inputs:
        probability - An Array of floats or a single float. If it is a single
        float, the same probability will be applied to every photon. If it is 
        an array, it must be the same the length as the Rays object, and the
        probability at index 0 will be applied to the first photon, repeated 
        for the rest of the array.
        
        Outputs:
        Nothing
        
        Notes:
        -Generates probabilities using the fact that np.random.rand() has a 20%
        chance to create a number below .2
        '''
        try:
            l = len(probability)
        except:
            probability = np.array([probability] * len(self))
        if len(probability) != len(self):
            raise ValueError('Input Array must be the same length as Rays Object')
        
        # With weighting, probability remove is considered differently, no
        # photons are actually removed
        if considerweights:
            if self.weighting:
                self.weights *= probability
                return
            
        choice = np.random.rand(len(self))
        needToRemove = choice > probability
        self.remove(needToRemove)
    
    
    ## Tag Functions:
    # These functions deal with 'Tags', which are truth arrays which are used
    # to mark certain photons in a Rays object.
    # Tags take the form of tuples: (tagname,trutharray)
    
    def addTag(self,tagname,trutharray):
        '''
        Function addTag:
        Creates a new tag for the Rays object
        
        Inputs:
        tagname - The name that the tag will be referenced with in the future
        trutharray - The array that tags individual photons in the Rays object
        
        Outputs:
        Nothing
        
        Notes:
        -Nominally, the trutharray is a Numpy array of Booleans, but an integer
        array will suffice where 0s are False and everything else is True
        -The input trutharray must be the same length as the Rays object
        -Trying to use a tagname that already exists will raise an error,
        the existing tag must be removed before a new one can be added
        -The tagname cannot start with '!' or '~', this will cause issues 
        in the combineTags method
        '''
        arr = self.getTag(tagname)
        if arr is not None:
            raise Exception('That tag already exists, remove the existing tag before adding a new one')
        
        if (tagname[0] == '!' or tagname[0] == '~'):
            raise ValueError('The first character of a Tag Name cannot be ! or ~')
        if (len(trutharray) == self.__len__()):
            self.tags.append((tagname,trutharray))
        else:
            raise Exception('Tag Array must be the same length as this Rays object')
    
    
    def removeTag(self,tagname):
        '''
        Function removeTag:
        Removes a tag for the Rays object
        
        Inputs:
        tagname - The name that the tag that you want to remove
        
        Outputs:
        Nothing
        
        Notes:
        -Raises an error if the tagname does not exist for this Rays object
        '''
        for i in range(len(self.tags)):
            if (tagname in self.tags[i]):
                del self.tags[i]
                return
        raise ValueError('That tag does not exist in this Rays object')
    
    
    def printTags(self):
        '''
        Function printTags:
        Prints the names of all Tags for this Rays object
        
        Inputs:
        Nothing
        
        Outputs:
        Nothing
        
        Notes:
        -The name of a Tag object is stored as the first element in its Tuple
        '''
        for i in self.tags:
            print(i[0])
    
    
    def getTag(self, tagname):
        '''
        Function getTag:
        Returns the trutharray for a Tag in this Rays object
        
        Inputs:
        tagname - The name of the tag that the user whishes to access
        
        Outputs:
        trutharray - The trutharray associated with that tagname
        
        Notes:
        -Returns None if the tag does not exist in this Rays object
        '''
        for i in range(len(self.tags)):
            if (tagname in self.tags[i]):
                return self.tags[i][1]
        return None
    
    
    def combineTags(self,tags,delim=None,orcombination=True):
        '''
        Function combineTags:
        Combines the trutharray of several tags so that multiple tags can be
        passed to a single function
        
        Inputs:
        tags - A string or a list of strings. If it is a list, each element 
        should be the name of a tag. If it is a string, the string should
        contain each tagname separated by a specified delimiting character. 
        If there is only one tagname in the string, leaving delim as None 
        is fine, but if there is more than one tag, delim should be specified.
        delim - The delimiting character that separates tagnames in a string
        orcombination - The method of combining tags, if True, the boolean 
        arrays will be combined using np.logical_or, if False, the boolean 
        arrays will be combined using np.logical_and
        
        Outputs:
        outputarray - The combined trutharray
        
        Notes:
        - The opposite (logical not) of an array can be used by this function.
        To do so, use the '!' or '~' characters before the tagname
        '''
        # delim cannot be ~ or !, this would interfere with opposite tagarrays
        if (delim == '~' or delim == '!'):
            raise ValueError("Delimiter cannot be '~' or '!'")
        
        # Start the output array as False, because other arrays will be added
        # using or statements
        outputarray = np.zeros(self.__len__())
        # Split the string based on delimiter if it is a string
        # Splitting based on None makes a singular-element array containing the
        # original string
        if type(tags) == type('this is a string'):
            tags = tags.split(delim)
            
        for tag in tags:
            # This block deals with tags that are specified as opposite
            if (tag[0] == '!' or tag[0] == '~'):
                trutharray = self.getTag(tag[1:])
                # Trutharray will be None if the tag is not in the Rays object
                if trutharray is None:
                    continue
                else:
                    # Note that the logical_not is passed into the outputarray
                    if orcombination:
                        outputarray = np.logical_or(outputarray,np.logical_not(trutharray))
                    else:
                        outputarray = np.logical_and(outputarray,np.logical_not(trutharray))
            # This block deals with tag not specified as opposite
            else:
                trutharray = self.getTag(tag)
                if trutharray is None:
                    continue
                else:
                    if orcombination:
                        outputarray = np.logical_or(outputarray,trutharray)
                    else:
                        outputarray = np.logical_and(outputarray,trutharray)
        
        return outputarray
    
    
    ## Parameter Functions:
    # Parameters are similar to tags but don't function as truth arrays.
    # Each element in a parameter array should be a float or integer.
    # A paramter can be passed to a function when that function requires 
    # imformation from each photon.
    # e.g: Some functions require information about the energy of each photon
    
    def addParam(self,paramname,paramarray):
        '''
        Function addParam:
        Creates a new parameter for the Rays object
        
        Inputs:
        tagname - The name that the parameter will be referenced with in the future
        trutharray - The array of values of that parameter for each photon in
        the Rays object
        
        Outputs:
        Nothing
        
        Notes:
        -The input paramarray must be the same length as the Rays object
        '''
        arr = self.getParam(paramname)
        if arr is not None:
            raise Exception('That Parameter already exists, remove the existing Parameter before adding a new one')
            
        if (len(paramarray) == self.__len__()):
            self.params.append((paramname,paramarray))
        else:
            raise Exception('Parameter Array must be the same length as this Rays object')
    
    def removeParam(self,paramname):
        '''
        Function removeParam:
        Removes a parameter from the Rays object
        
        Inputs:
        paramname - The name that the parameter that you want to remove
        
        Outputs:
        Nothing
        
        Notes:
        -Raises an error if the Parameter does not exist for this Rays object
        '''
        for i in range(len(self.params)):
            if (paramname in self.params[i]):
                del self.params[i]
                return
        raise ValueError('That Parameter does not exist in this Rays object')
    
    
    def printParams(self):
        '''
        Function printParam:
        Prints the names of all parameters for this Rays object
        
        Inputs:
        Nothing
        
        Outputs:
        Nothing
        
        Notes:
        -The name of a Parameter object is stored as the first element in its Tuple
        '''
        for i in self.params:
            print(i[0])
    
    
    def getParam(self, paramname):
        '''
        Function getTag:
        Returns the paramarray for a Parameter in this Rays object
        
        Inputs:
        paramname - The name of the parameter that the user whishes to access
        
        Outputs:
        valuearray - The array that stores values for this specific parameter
        
        Notes:
        -Returns None if the parameter does not exist in this Rays object
        '''
        for i in range(len(self.params)):
            if (paramname in self.params[i]):
                return self.params[i][1]
        return None
    
    
    ## Addition Functions:
    # Control what happens when you try to add two Rays objects
    
    def __add__(self,other):
        '''
        Function __add__:
        Magic method overrides the + operator for Ray objects
        
        Inputs:
        self - The first Rays object in the expression
        other - The second Rays object in the expression
        
        Outputs:
        newrays - The sum of the two input Rays objects
        
        Notes:
        -The x,y,z,l,m,n,ux,uy,uz arrays from the two Rays objects are 
        concatenated to form the new Rays object
        -If the same tag or parameter exists in both Rays objects, their
        respective arrays are concatenated for the output Rays object
        -If Ray object A (one of the two inputs) has a tag or parameter that
        Ray object B (the other input) does not have, a default array of
        length len(B) is created and appended to A's tag
        -The default arrays are entirely zeros (representing False for tag
        arrays), if you wish to use different default values, use the 
        alternate add method (not __add__)
        -Called with the synax:
        >> C=A+B
        where A and B are Rays objects
        '''
        return self.add(other)
    
    
    def __iadd__(self,other):
        '''
        A wrapper for the __add__ method, has the same functionality,
        but is called with the syntax:
        >> A+=B
        where A and B are Rays objects, and A becomes the output of __add__
        '''
        return self.__add__(other)
    
    def add(self,other,defaulttag=0,defaultparam=0):
        '''
        Function __add__:
        Magic method overrides the + operator for Ray objects
        
        Inputs:
        self - The first Rays object in the expression
        other - The second Rays object in the expression
        
        Outputs:
        newrays - The sum of the two input Rays objects
        
        Notes:
        -The x,y,z,l,m,n,ux,uy,uz arrays from the two Rays objects are 
        concatenated to form the new Rays object
        -If the same tag or parameter exists in both Rays objects, their
        respective arrays are concatenated for the output Rays object
        -If Ray object A (one of the two inputs) has a tag or parameter that
        Ray object B (the other input) does not have, a default array of
        length len(B) is created and appended to A's tag
        -This method allows you to choose the contents of your default arrays
        -The default arguments are 0 for both tags and parameters, but these
        can be specified using the defaulttag and defaultparameter arguments,
        respectively
        -defaulttag and defaultparam arguments must be numeric (0 will be False
        for the tag arrays)
        -If one of the two Rays objects has weighting, it will be added to the
            other Rays object
        -If both Rays objects have weighting, their weights will be concatenated
        -If neither Rays objects have weighting, the resultant Rays object will
            also not have weighting
        '''
        tagstr = str(defaulttag)
        paramstr = str(defaultparam)
        if (not tagstr.isnumeric() or not paramstr.isnumeric()):
            raise ValueError('Defaulttag and DefaultParam must both be numeric')
        
        selflen = len(self)
        otherlen = len(other)
        
        # Concatenate basic arrays
        xarr = np.concatenate((self.x,other.x))
        yarr = np.concatenate((self.y,other.y))
        zarr = np.concatenate((self.z,other.z))
        larr = np.concatenate((self.l,other.l))
        marr = np.concatenate((self.m,other.m))
        narr = np.concatenate((self.n,other.n))
        uxarr = np.concatenate((self.ux,other.ux))
        uyarr = np.concatenate((self.uy,other.uy))
        uzarr = np.concatenate((self.uz,other.uz))
        
        # Create new Rays object with new arrays
        newrays = Rays(x=xarr,y=yarr,z=zarr,l=larr,m=marr,n=narr,ux=uxarr,uy=uyarr,uz=uzarr)
        
        # Iterate through self's tags
        for tag in self.tags:
            # Check if other has the same tag
            othertag = other.getTag(tag[0])
            if othertag is None: 
                # If it doesn't have the same tag, create an array of value
                # defaulttag
                othertag = np.ones(otherlen) * defaulttag
            # Add the new tag to the new Rays object
            newrays.addTag(tag[0],np.concatenate((self.getTag(tag[0]),othertag)))
        
        # Iterate through other's tags
        for tag in other.tags:
            # Check if the tag already exists in newrays
            currenttag = newrays.getTag(tag[0])
            if currenttag is not None:
                # If the tag alreade exists in newrays, skip this tag
                continue
            # Check if the current tag exists in self
            selftag = self.getTag(tag[0])
            if selftag is None:
                # If the current tag does not exist in self, create an array for
                # it using defaulttag
                selftag = np.ones(selflen) * defaulttag
            # Add the new tag to the Rays object
            newrays.addTag(tag[0],np.concatenate((selftag,other.getTag(tag[0]))))
        
        # Does the same for parameters that it just did for tags
        for param in self.params:
            otherparam = other.getParam(param[0])
            if otherparam is None: 
                otherparam = np.ones(otherlen) * defaultparam
            newrays.addParam(param[0],np.concatenate((self.getParam(param[0]),otherparam)))
        
        for param in other.params:
            currentparam = newrays.getParam(param[0])
            if currentparam is not None:
                continue
            selfparam = self.getParam(param[0])
            if selfparam is None:
                selfparam = np.ones(selflen) * defaultparam
            newrays.addParam(param[0],np.concatenate((selfparam,other.getParam(param[0]))))
        
        # Check to see if either rays objects has weighting
        if self.weighting and not other.weighting:
            otherweights = np.ones(len(other))
            newweights = np.concatenate((self.weights,other.weights))
        elif not self.weighting and other.weighting:
            selfweights = np.ones(len(self))
            newweights = np.concatenate((selfweights,other.weights))
        elif self.weighting and other.weighting:
            newweights = np.concatenate((self.weights,other.weights))
        else:
            newweights = None
        if newweights is not None:
            newrays.weighting = True
            newrays.weights = newweights
        
        # If both objects have wavelengths
        if (self.wave is not None and other.wave is not None):
            newrays.wave = np.concatenate((self.wave,other.wave))
        #TODO: Currently this supposes that the only way wave could be None is for an empty Rays object, improve in the future
        # If only one has wavelengths
        elif (self.wave is not None and other.wave is None):
            newrays.wave = self.wave
        elif (self.wave is None and other.wave is not None):
            newrays.wave = other.wave
        # If both have wavelengths, do nothing
        
        # If both objects have order
        if (self.order is not None and other.order is not None):
            newrays.order = np.concatenate((self.order,other.order))
        #TODO: Currently this supposes that the only way order could be None is for an empty Rays object, improve in the future
        # If only one has order
        elif (self.order is not None and other.order is None):
            newrays.order = self.order
        elif (self.order is None and other.order is not None):
            newrays.order = other.order
        # If both have order, do nothing
        
        
        
        
        return newrays
    
    
    ## Efficiency Calculations:
    # Functions needed to calculate the efficiency of optical systems
    # Used to keep track the number of photons in a Rays object
    
    def starteffcalc(self):
        '''
        Function starteffcalc:
        Saves the current length of the array, used before photons will be
        removed to serve as a starting point in an efficiency calculation
        
        Inputs:
        Nothing
        
        Outputs:
        Nothing
        
        Notes:
        -Calling this function again will overwrite the old starting point
        '''
        self.startingeff = self.__len__()
    
    def checkeffcalc(self):
        '''
        Function checkeffcalc:
        Returns the paramarray for a Parameter in this Rays object
        
        Inputs:
        Nothing
        
        Outputs:
        startingeff - The saved photon length from a previous call of 
        starteffcalc()
        
        Notes:
        -Returns None starteffcalc() has not yet been called for this Rays
        object
        '''
        return self.startingeff
    
    def geteff(self,display=False,header='Efficiency Calculation'):
        '''
        Function geteff:
        Performs an efficiency calculation
        
        Inputs:
        display - This argument should be true if the user wants the function
        to print out information about the efficiency calculation
        header - A user defined header to differentiate several efficiency
        calculations, if None is specified, the header will be 'Efficiency
        Calculation'
        
        Outputs:
        efficiency - The fraction of photons that remain since the last call
        of starteffcalc()
        
        Notes:
        -Calling this function doesn't update self.startingeff, so efficiency
        calculations can be performed at several ending points from the same
        starting point
        -An example of the display (using default header) is:
        Efficiency Calculation
        ----------------------
        8 photons started
        8 photons remain
        Efficiency: 1.0
        ======================
        '''
        if self.startingeff is None:
            # Raise an exception if the efficiency calculation has not been started
            raise Exception('Efficiency Calculation needs to be started with starteffcalc()')
        starteff = self.startingeff
        endeff = self.__len__()
        if display:
            print(header)
            print('----------------------')
            print(str(starteff) + ' photons started')
            print(str(endeff) + ' photons remain')
            print('Efficiency: ' + str(endeff / starteff))
            print('======================')
        return endeff / starteff
    
    

    
    ## Graphing Functions:
    
    
    def histogram(self,param='x',     
                  bins=None,range=None,density=None,weights=None,
                  cumulative=False,bottom=None,histtype='bar',align='mid',
                  orientation='vertical',rwidth=None,log=False,color=None,
                  label=None,stacked=False,normed=None):
        
        if param.lower() == 'x':
            arg = self.x
            label = 'X'
        elif param.lower() == 'y':
            arg = self.y
            label = 'Y'
        elif param.lower() == 'z':
            arg = self.z
            label = 'Z'
        elif param.lower() == 'l':
            arg = self.l
            label = 'L'
        elif param.lower() == 'm':
            arg = self.m
            label = 'M'
        elif param.lower() == 'n':
            arg = self.n
            label = 'N'
        elif param.lower() == 'ux':
            arg = self.ux
            label = 'UX'
        elif param.lower() == 'uy':
            arg = self.uy
            label = 'UY'
        elif param.lower() == 'uz':
            arg = self.uz
            label = 'UZ'
        elif param.lower() == 'pos' or param.lower() == 'position':
            arg = np.sqrt(self.x**2 + self.y**2 + self.z**2)
            label = 'Position Magnitude'
        elif param.lower() == 'dir' or param.lower() == 'direction':
            arg = np.sqrt(self.l**2 + self.m**2 + self.n**2)
            label = 'Direction Magnitude'
        elif param.lower() == 'normal' or param.lower() == 'surfacenormal':
            arg = np.sqrt(self.ux**2 + self.uy**2 + self.uz**2)
            label = 'Surface Normal Magnitude'
        else:
            raise ValueError('param argument is not a recognized Parameter for the histogram')
        
        plt.figure()
        plt.hist(arg,bins=bins,range=range,density=density,weights=weights,
                     cumulative=cumulative,bottom=bottom,histtype=histtype,
                     align=align,orientation=orientation,rwidth=rwidth,
                     log=log,color=color,label=label,stacked=stacked,
                     normed=normed)
        plt.xlabel(label)
        plt.show()
    
    
    def scatter2d(self,horiz='x',vert='y',s=None,c=None,marker=None,cmap=None,
                       norm=None,vmin=None,vmax=None,alpha=None,linewidths=None, 
                       verts=None,edgecolors=None):
        
        if horiz.lower() == 'x':
            harg = self.x
            hlabel = 'X'
        elif horiz.lower() == 'y':
            harg = self.y
            hlabel = 'Y'
        elif horiz.lower() == 'z':
            harg = self.z
            hlabel = 'Z'
        elif horiz.lower() == 'l':
            harg = self.l
            hlabel = 'L'
        elif horiz.lower() == 'm':
            harg = self.m
            hlabel = 'M'
        elif horiz.lower() == 'n':
            harg = self.n
            hlabel = 'N'
        elif horiz.lower() == 'ux':
            harg = self.ux
            hlabel = 'UX'
        elif horiz.lower() == 'uy':
            harg = self.uy
            hlabel = 'UY'
        elif horiz.lower() == 'uz':
            harg = self.uz
            hlabel = 'UZ'
        elif horiz.lower() == 'pos' or horiz.lower() == 'position':
            harg = np.sqrt(self.x**2 + self.y**2 + self.z**2)
            hlabel = 'Position Magnitude'
        elif horiz.lower() == 'dir' or horiz.lower() == 'direction':
            harg = np.sqrt(self.l**2 + self.m**2 + self.n**2)
            hlabel = 'Direction Magnitude'
        elif horiz.lower() == 'normal' or horiz.lower() == 'surfacenormal':
            harg = np.sqrt(self.ux**2 + self.uy**2 + self.uz**2)
            hlabel = 'Surface Normal Magnitude'
        else:
            raise ValueError('horiz argument is not a recognized Parameter for the histogram')
            
        if vert.lower() == 'x':
            varg = self.x
            vlabel = 'X'
        elif vert.lower() == 'y':
            varg = self.y
            vlabel = 'Y'
        elif vert.lower() == 'z':
            varg = self.z
            vlabel = 'Z'
        elif vert.lower() == 'l':
            varg = self.l
            vlabel = 'L'
        elif vert.lower() == 'm':
            varg = self.m
            vlabel = 'M'
        elif vert.lower() == 'n':
            varg = self.n
            vlabel = 'N'
        elif vert.lower() == 'ux':
            varg = self.ux
            vlabel = 'UX'
        elif vert.lower() == 'uy':
            varg = self.uy
            vlabel = 'UY'
        elif vert.lower() == 'uz':
            varg = self.uz
            vlabel = 'UZ'
        elif vert.lower() == 'pos' or vert.lower() == 'position':
            varg = np.sqrt(self.x**2 + self.y**2 + self.z**2)
            vlabel = 'Position Magnitude'
        elif vert.lower() == 'dir' or vert.lower() == 'direction':
            varg = np.sqrt(self.l**2 + self.m**2 + self.n**2)
            vlabel = 'Direction Magnitude'
        elif vert.lower() == 'normal' or vert.lower() == 'surfacenormal':
            varg = np.sqrt(self.ux**2 + self.uy**2 + self.uz**2)
            vlabel = 'Surface Normal Magnitude'
        else:
            raise ValueError('vert argument is not a recognized Parameter for the histogram')
        
        plt.figure()
        plt.scatter(harg,varg,s=s,c=c,marker=marker,cmap=cmap,norm=norm,
                    vmin=vmin,vmax=vmax,alpha=alpha,linewidths=linewidths, 
                    verts=verts,edgecolors=edgecolors)
        plt.xlabel(hlabel)
        plt.ylabel(vlabel)
        plt.show()
        
        
    def scatter3d(self,type='pos',s=None,c=None,marker=None,cmap=None,
                       norm=None,vmin=None,vmax=None,alpha=None,linewidths=None, 
                       verts=None,edgecolors=None):
        # Import statement not at header because I don't want this library
        # to be a requirement for the whole program
        try:
            from mpl_toolkits.mplot3d import Axes3D
        except:
            raise Exception('The scatter3d method requires the mplot3d API')
        
        if type.lower() == 'pos' or type.lower() == 'position':
            xarg = self.x
            yarg = self.y
            zarg = self.z
            xlabel = 'X'
            ylabel = 'Y'
            zlabel = 'Z'
        elif type.lower() == 'dir' or type.lower() == 'direction':
            xarg = self.l
            yarg = self.m
            zarg = self.n
            xlabel = 'L'
            ylabel = 'M'
            zlabel = 'N'
        elif type.lower() == 'normal' or type.lower() == 'surfacenormal':
            xarg = self.ux
            yarg = self.uy
            zarg = self.uz
            xlabel = 'UX'
            ylabel = 'UY'
            zlabel = 'UZ'
        else:
            raise ValueError('type argument is not a recognized Parameter for the histogram')
        
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        p=ax.scatter(xarg,yarg,zarg,c=c,marker=marker,cmap=cmap,norm=norm,
                   vmin=vmin,vmax=vmax,alpha=alpha,linewidths=linewidths, 
                   verts=verts,edgecolors=edgecolors)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        if c is not None:
            plt.colorbar(p)
        plt.show()
    