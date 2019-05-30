import numpy as np
# from prtp.Rays import Rays

class Source:
    
    def __init__(self, num = 0, wave = None, order = None):
        self.type = type
        self.num = num
        self.wave = wave
        self.order = order
        self.rays = None
    
    def addMisalignment(self, dx = 0,dy = 0, dz = 0, dl = 0, dm = 0, dn = 0):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.dl = dl
        self.dm = dm
        self.dn = dn
    
    def misalign(self,rays):
        try:
            # Add misalignments if they've been defined
            rays.translate(self.dx,self.dy,self.dz)
            rays.rotate(self.dl,self.dm,self.dn)
            rays.rotatenormal(self.ux,self.uy,self.uz)
        except:
            # This executes if misalignments have not been defined
            pass
    
    def addParams(self,rays):
        if (self.wave is not None):
            if type(self.wave) == np.ndarray:
                rays.addParam('Wavelength',self.wave)
            else:
                rays.addParam('Wavelength',np.array([self.wave]*len(rays)))
        if (self.order is not None):
            if type(self.wave) == np.ndarray:
                rays.addParam('Order',self.order)
            else:
                rays.addParam('Order',np.array([self.order]*len(rays)))
    
    def getRays():
        if self.rays is None:
            self.generateRays()
        return self.rays