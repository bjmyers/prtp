import numpy as np
import astropy.units as u
from prtp.Sources import Subannulus
from prtp.WolterOptic import WolterModule
from prtp.Rays import Rays
from prtp.Instrument import Instrument

# Define the inner and outer radii of the subannulus
rp_front = 195. * u.mm
rp_back = 165. * u.mm

# Generate the subannulus source
r = Subannulus(1000,rp_back, rp_front, np.radians(30.)*u.rad,wave=0.83401*u.nm,order=0)

# Generate the list of r0s for the various mirrors, these values taken from OGRE
r0s = np.array([165., 167.5503, 170.1193, 172.7023,
                  175.3143, 177.9404, 180.5859, 183.2509,
                  185.9355, 188.6398, 191.3640, 194.1083]) * u.mm
# NOTE: There are 12 mirrors in this module, so all of the other specifications
# need to be arrays of length 12

# Generate lists of z0, axial_lengths, and mirror_separations, values taken from OGRE
# Since every mirror will have the same value for these parameters, we can generate
# the arrays using this syntax:
z0s = np.ones(12) * 3500. * u.mm
ax_lens = np.ones(12) * 100 * u.mm
mir_seps = np.ones(12) * 5 * u.mm

# Generate the Wolter Module
wm = WolterModule(r0=r0s,z0=z0s,axial_length=ax_lens,mirror_sep=mir_seps)
# NOTE: Since we did not define a position (x,y,z) or direction vector (nx,ny,nz)
# This Wolter Module will use the default values:
# position = (0,0,0)
# direction = (0,0,1)

# Initialize the Instrument with the source
i = Instrument(r)
# Add the only component to the Instrument
i.addComponent(wm)
# Simulate the Rays through the instrument
i.simulate()
# Print efficiency information for the user
i.displayEfficiency()

# Access the final rays and plot them in 3D
rays = i.getRays()
rays.scatter3d()
# NOTE: The first plot shown will not show the individual mirrors obviously.
# Drag the plot until you're looking down the z-axis. Then you'll see the
# distinctive shapes of the concentric mirrors.