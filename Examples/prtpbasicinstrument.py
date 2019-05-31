import numpy as np
import prtp
import prtp.conicsolve
from prtp.Rays import Rays
from prtp.WolterOptic import WolterTypeOne
from prtp.Grating import Grating
from prtp.Instrument import Instrument
from prtp.Modification import Modification
from prtp.Sources import Subannulus
import astropy.units as u

# Define Wolter-I parameters.
r0 = (165. * u.mm).to('parsec')  # Radius at the intersection plane of primary/secondary.
z0 = 3500. * u.mm 
mirror_len = 100. * u.mm
mirror_sep = 5. * u.mm

d = 160. * u.nm  # Groove period @ 3300 mm [nm]
L = 3250. * u.mm # Center of grating [mm]
d *= L  / (3300 * u.mm)  # Find groove period at center of grating.

# Define inner & outer radii to create rays.
rp_front = prtp.conicsolve.primrad(z0 + mirror_sep/2 + mirror_len, r0, z0)
rp_back = prtp.conicsolve.primrad(z0 + mirror_sep/2, r0, z0)

# Define initial rays in subannulus.
source = Subannulus(1000,rp_back, rp_front, np.radians(30.)*u.rad,wave=0.83401*u.nm,order=0)

# Define Rotation Modification to change the Rays so that diffraction
# occurs in the x-axis
def rot(rays,cw):
    rays.transform(0, 0, 0, 0, 0, np.radians(90.))
rotation = Modification(rot)

# Define Wolter Optic
wolter = WolterTypeOne(r0=r0,z0=z0)

# Define Grating (values from old code)
grat = Grating(0.*u.mm,151.86466758*u.mm,3247.56521956*u.mm,
            0.,0.99978888,-0.0205472,
            0.01518378,0.02054483,0.99967363,
            l=100*u.mm,w=100*u.mm,d=d,radial=True,fdist=3250*u.mm)

def func(rays,cw):
    rays.beckmann_scatter(0,0,1.48e-5)
scatter = Modification(func)

i = Instrument(source)
i.addComponent(rotation)
i.addComponent(wolter)
i.addComponent(scatter)
i.addComponent(grat)
i.simulate()
i.displayEfficiency()

rays = i.getRays()

rays.focusX()
rays.scatter3d()


#PyXFocus: 1M photons thru basic instrument in 3.89868 s
#New PRTP: 1M photons thru basic instrument in 5.03021 s