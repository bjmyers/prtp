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
r0 = (165. * u.mm)  # Radius at the intersection plane of primary/secondary.
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
source = Subannulus(10000,rp_back, rp_front, np.radians(30.)*u.rad,wave=0.83401*u.nm,order=0)

# Define Wolter Optic
# Make sure to add Beckmann Scattering, will be added after the Primary Mirror
wolter = WolterTypeOne(r0=r0,z0=z0,beckmann_scatter=True,ripple=1.48e-5)

# Define Grating (values from old code)
grat = Grating(0.*u.mm,151.86466758*u.mm,3247.56521956*u.mm,
            0.,0.99978888,-0.0205472,
            -0.01518378,-0.02054483,-0.99967363,
            l=100*u.mm,w=100*u.mm,d=d,radial=True,fdist=3250*u.mm)

# Initialize the instrument and add components
i = Instrument(source)
i.addComponent(wolter)
i.addComponent(grat)
i.addFocus()

# Simulate the Rays through the instrument
i.simulate()

# Print efficiency information for the user
i.displayEfficiency()

# Access the final rays
rays = i.getRays()

# Plot the Rays
rays.scatter2d()
