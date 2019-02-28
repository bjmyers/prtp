import numpy as np
import matplotlib.pyplot as plt
import sources
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

class Rays:
    
    ## Creation Functions:
    
    
    def __init__(self,x=[],y=[],z=[],l=[],m=[],n=[],ux=[],uy=[],uz=[]):
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.l = np.array(l)
        self.m = np.array(m)
        self.n = np.array(n)
        self.ux = np.array(ux)
        self.uy = np.array(uy)
        self.uz = np.array(uz)
        # This call to the __len__ method makes sure the inputs were all of the
        # same length
        self.__len__()
        self.tags = []
        self.params = []
    
    
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
                                           zset,rin,rout,tmin,tmax,num,lscat)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
    
    @classmethod
    def rectbeam(cls,xhalfwidth,yhalfwidth,num):
        opd,x,y,z,l,m,n,ux,uy,uz = sources.rectbeam(
                                           xhalfwidth,yhalfwidth,num)
        return cls(x,y,z,l,m,n,ux,uy,uz)
    
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
    
    
    def combineTags(self,tags,delim=None):
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
        
        Outputs:
        outputarray - The combined trutharray
        
        Notes:
        - The opposite (logical not) of an array can be used by this function.
        To do so, use the '!' or '~' characters before the tagname
        '''
        # delim cannot be ~ or !, this would interfere with opposite tagarrays
        if (delim == '!' or delim == '!'):
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
                    outputarray = np.logical_or(outputarray,np.logical_not(trutharray))
            # This block deals with tag not specified as opposite
            else:
                trutharray = self.getTag(tag)
                if trutharray is None:
                    continue
                else:
                    outputarray = np.logical_or(outputarray,trutharray)
        
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
        
        return newrays
    
    
    ## Basic Motion Functions:
    
    def translate(self,dx=0,dy=0,dz=0):
        self.x += dx
        self.y += dy
        self.z += dz
    
    
    def rotate(self,dl=0,dm=0,dn=0):
        self.l += dl
        self.m += dm
        self.n += dn
    
    
    def rotatenormal(self,dux=0,duy=0,duz=0):
        self.ux += dux
        self.uy += duy
        self.uz += duz
    
    
    def move(self,x=None,y=None,z=None,
                  l=None,m=None,n=None,
                  ux=None,uy=None,uz=None):
        size = len(self.x)
        if x is not None:
            self.x = np.ones(size) * x
        if y is not None:
            self.y = np.ones(size) * y
        if z is not None:
            self.z = np.ones(size) * z
        if l is not None:
            self.l = np.ones(size) * l
        if m is not None:
            self.m = np.ones(size) * m
        if n is not None:
            self.n = np.ones(size) * n
        if ux is not None:
            self.ux = np.ones(size) * ux
        if uy is not None:
            self.uy = np.ones(size) * uy
        if uz is not None:
            self.uz = np.ones(size) * uz
    
    
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
        ax.scatter(xarg,yarg,zarg,c=c,marker=marker,cmap=cmap,norm=norm,
                   vmin=vmin,vmax=vmax,alpha=alpha,linewidths=linewidths, 
                   verts=verts,edgecolors=edgecolors)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        plt.show()
    
    
x = Rays.circularbeam(5,3)
x.addParam('tag1',[0,1,0])
x.addParam('sharedTag',[1,1,0])
y = Rays.circularbeam(3,5)
y.addParam('tag2',[1,1,1,1,1])
y.addParam('sharedTag',[0,1,1,1,1])
    
    
    
    
    
    
    
    
    
    
    
    
    
    