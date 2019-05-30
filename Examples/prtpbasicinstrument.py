import numpy as np
import prtp
import prtp.conicsolve
from prtp.Rays import Rays
from prtp.WolterTypeOne import WolterTypeOne
from prtp.Grating import Grating
from prtp.Instrument import Instrument
from prtp.Modification import Modification

# Define Wolter-I parameters.
r0 = 165.  # Radius at the intersection plane of primary/secondary.
z0 = 3500. # 
mirror_len = 100.
mirror_sep = 5.

d = 160.  # Groove period @ 3300 mm [nm]
L = 3250.  # Center of grating [mm]
d *= L  / 3300  # Find groove period at center of grating.

# Define inner & outer radii to create rays.
rp_front = prtp.conicsolve.primrad(z0 + mirror_sep/2 + mirror_len, r0, z0)
rp_back = prtp.conicsolve.primrad(z0 + mirror_sep/2, r0, z0)

# Define initial rays in subannulus.
rays = Rays.subannulus(rp_back, rp_front, np.radians(30.), 1000)

# Rotate 90 degrees so that diffraction occurs in x-axis.
rays.transform(0, 0, 0, 0, 0, np.radians(90.))

w = WolterTypeOne(r0=r0,z0=z0)
g = Grating(0.,151.86466758,3247.56521956,
            0.,0.99978888,-0.0205472,
            0.01518378,0.02054483,0.99967363,
            l=100,w=100,d=d,radial=True,fdist=3250)

def func(rays,cw):
    rays.beckmann_scatter(0,0,1.48e-5)
    
m = Modification(func)

i = Instrument(rays,waves=0.83401,orders=-2,weighting=False)
i.addComponent(w)
i.addComponent(m)
i.addComponent(g)
i.simulate()
i.displayEfficiency()


rays.focusX()
rays.scatter3d()


#PyXFocus: 1M photons thru basic instrument in 3.89868 s
#New PRTP: 1M photons thru basic instrument in 5.03021 s