
.. _modification-top:

Modification
================================

A Modification is essentially a "catch-all" component which allows the user to modify rays within an Instrument in the middle of an Instrument.

For example, suppose the user had a special type of scattering the Rays were supposed to undergo after leaving a Wolter Module. This can not be done by default by any Components. So doing this would require the user to stop the Instrument after tracing to the Wolter Module, get the rays, add the scattering, package the rays back into a source, and feed them into another Instrument containing the rest of the components.

A Modification would have allowed the user to add the scattering without interrupting the Instrument's simulation.

Creating a Modification
-------------------------

Modifications are created with a single argument:

* function - A function that takes in Rays and modifies them in place, see examples for more information
   * The function takes in two arguments: rays and considerweights. Rays is the Rays object you should modify, and considerweights can be used if weighted rays should be treated differentley. Though this argument is not often used, it must be included in the function header. 
   * Any return statements will be ignored, the rays should be modified in place.

When rays are traced to a Modification using the trace() function, the function defined in initialization is called.

Example
------------

Suppose you wanted to simulate a special type of scattering where the l of each photon is offset by a Gaussian.

Before defining the Modification, we must defined the scattering function. To work within a Modification, the function must take in two arguments, rays and considerweights. However, considerweights will have no effect on our function.

.. code-block:: python

   import numpy as np

   def scatter(rays, considerweights):
   
      # Add the Gaussian with std. of .05
      rays.l += np.random.normal(0,.05,len(rays))

      # Re-normalize the direction length
      rays.n = np.sqrt(rays.l**2 + rays.m**2)

      # Do not return rays, we have modified them correctly


Now that we have our function defined, we can create a Modification:

.. code-block:: python

   from prtp.Modification import Modification

   m = Modification(scatter)

Now our Modification m can be added to an Instrument.






