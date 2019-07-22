
.. _instrument-top:

Instruments
===================

Instruments are objects that allow the user to simulate several components in sequence. Instrument objects allow for the reporting of detailed efficiency information and automated misalignment tests.

Instruments are unlike Combinations. Combinations contain a list of components but will trace the rays to all of the components at once, so that a single photon will only ever be traced to one of the components. Instruments, however, trace rays to the components in sequence, so the photons will be traced to the first component, then the second, and so on, such that the final surviving photons have been successfully traced to every component in the Instrument.

Creating an Instrument
-------------------------

Instrument objects take the following arguments:

* source - A Source object which will produce the starting rays for this Instrument.
   * If you wish to send rays from one Instrument directly into the next. You can do so by creating the first instrument, simulating it, then passing that instrument as the source argument for the second instrument. See the example at the end of this seciton.
* considerweights - A boolean that determines if the photons should be weighted or not for this simulation. If True, the photons will be weighted.

Example 1: A basic Instrument Definition

.. code-block:: python

   from prtp.Sources import CircularBeam
   from prtp.Instrument import Instrument
   import astropy.units as u

   s = CircularBeam(rad=2*u.mm)

   inst = Instrument(s)


Example 2: Transferring rays from one Instrument to another

.. code-block:: python

   from prtp.Sources import CircularBeam
   from prtp.Instrument import Instrument
   import astropy.units as u

   s = CircularBeam(rad=2*u.mm)

   inst = Instrument(s)

   # Add components here

   inst.simulate()

   # Create a new instrument that will start with the rays produced
   # by the previous instrument
   nextinst = Instrument(inst)

   # Add components to second instrument here

   nextinst.simulate()

.. warning::

   Before defining an instrument with a previous instrument. The previous instrument must be simulated. Otherwise there will be no output rays to transfer and an error will be raised.

:ref:`Back to Top<instrument-top>`

Adding Components
-------------------

Components are added to an Instrument in the same way that they are added to a Combination. That is, they are added to the componentlist parameter with the function addComponent(). This function takes two arguments:

* comp - The Component you wish to add to the Instrument.
* index - The index at which you wish to add the Component. This argument defaults to None, in which case the Component will be appended to the end of the list. For example, if index was set to 3, the component will either be added at index 3 (with later components pushed back) or at the end of the componentlist if it had a length less than 3.

Example:

.. code-block:: python

   from prtp.Combination import Combination
   from prtp.CollimatorPlate import CollimatorPlate
   from prtp.Grating import Grating

   comb = Combination()

   c0 = CollimatorPlate()
   c1 = CollimatorPlate()
   c2 = CollimatorPlate()

   g = Grating()

   comb.addComponent(c0)

   # Without an index argument, the component will be added at the end
   comb.addComponent(c1)

   # If the index argument is too large, the component 
   # will be added at the end
   comb.addComponent(c2,index=100)

   # Add the Grating to be the second item in the list
   comb.addComponent(g,index=1)

:ref:`Back to Top<instrument-top>`


Adding a Focus
***************

Often, you will want to focus your rays after sending them through an instrument. Normally, this would be done by accessing the final rays and then calling focusX() on them. However, focusing can be performed within the Instrument using the addFocus() function. This function takes one argument:

* index - The index at which you want to add the focus. This argument follows the same rules as it does in the addComponent() function.

When this function is called, it adds a :ref:`Modification <modification-top>` object which takes the rays and calls the focusX() function on them.

Example:

.. code-block:: python

   from prtp.Sources import Circularbeam
   from prtp.Grating import Grating
   from prtp.Instrument import Instrument

   s = Circularbeam()

   g = Grating()

   inst = Instrument(s)
   inst.addComponent(g)
   inst.addFocus()

   inst.simulate()


Removing Components
----------------------

Components can be removed as easily as they can be added. This is done with the function removeComponent(), which takes the following arguments:

* index - The index of the component that you would like to remove.

Note: since this function uses the pop() function on the componentlist, it also returns the component which was removed.

Example:

.. code-block:: python

   from prtp.Sources import Circularbeam
   from prtp.Grating import Grating
   from prtp.Instrument import Instrument

   s = Circularbeam()

   g = Grating()

   inst = Instrument(s)

   inst.addComponent(g)

   inst.removeComponent(0)

   print(inst.componentlist)

This code block shows that the component list is now empty with the following output:

.. code-block:: python

   >>[]


Simulating The Components
---------------------------

Once you have added the components you would like to call (into an Instrument i, for example), you only need to call i.simulate() to trace the rays to each component in the Instrument.

The call to i.simulate() will generate rays from the source, then take those rays and trace them to every component in the componentlist, keeping those photons which interact with each component.

Calling i.simulate() again will generate new rays from the source and send those rays through the componentlist, it will not change the rays which were simulated previously.


Accessing Simulated Rays
-----------------------------

Once the Instrument has been simulated, the resultant rays can be accessed with the getRays() function, which takes no arguments.

Note that if the instrument has not yet been simulated, getRays() will return None.

:ref:`Back to Top<instrument-top>`

Displaying Simulation Information
------------------------------------

Efficiency Information
************************

Every component in an Instrument produces efficiency information when it is traced. This information includes how many rays were input, and how many rays were successfully traced. This information is collected by the Instrument object and can be displayed with the displayEfficiency() command, which takes no arguments. When the instrument from prtpbasicinstrument.py (in the Examples folder) produces its efficiency information, it produces the following output:

Method:                        Local Percent:   Global Percent:

Missed Wolter Optic            000.00000%       000.00000%

Missed Grating                 050.50000%       050.50000%

Failed to Reflect off Grating  000.00000%       000.00000%

Total Throughput: 49.50000%


^This output is better formatted in a Python shell.

Note that each method has two percentages, a local and a global percent. The local percent describes what percentage of the rays which attempted to be traced to the component were removed. The global percent describes the percent of the total photons (created by the source) which were removed by this method. We can see that 50.5% of the photons which were traced to the Grating were removed. Since no other photons were removed in this simulation, the global percentage is also 50.5%.

At the bottom is a total throughput, this describes the percentage of photons which made it through the entire instrument.

Spectral Resolution
***********************

After an Instrument has been simulated, the spectral resolution of the instrument can be determined using the spectralResolution() command, which takes no arguments. If the instrument has not yet been simulated, an error will be raised.

:ref:`Back to Top<instrument-top>`

Misalignment Tests
-------------------

These tests allow the user to automatically perform misalignments on components and compare the results. They will systematically misalign a given component, simulate the instrument, then record the FWHM of the photons with respect to a given variable. It is highly recommended that a focus be added to the end of the Instrument. Otherwise the FWHM of the photons has little meaning (in x, y, or z, that is).

All misalignment tests return two arrays. The first one containing the misalignment values, the second one containing the FWHM values.

These misalignment tests are still in the early stages of development, and thus may not always perform as expected.

Single Translation Test
*************************

This test takes a component and translates it along a specified axis. The function singleTranslateTest() takes the following arguments:

* index - The index of the component you would like to misalign, defaults to 0.
* min, max - The minimum and maximum translational values you would like to test. An equal number of steps will be chosen between these values to be tested. Default to -1 mm and 1 mm, respectively.
   * min and max must be in units of length, see the section on :ref:`Astropy Units <units-top>`
* num - The number of different tests you would like to perform. The misalignment values will be equally spaced between the min and max values. Defaults to 10.
* dim - The dimension along which you would like to translate. This argument should have a value of 1, 2, or 3 for the x, y, or z axis, respectively. Defaults to 1.
* plot - A boolean that tells the function if you would like to results to be automatically plotted. If True, a plot will be generated showing the full-width at half-max in the desired dimension as a function of misalignment value. Defaults to True.
* param - A string that tells the function the variable that you would like to plot the FWHM of. 
   * param has potential values: 'x', 'y', 'z', 'l', 'm', 'n', 'ux', 'uy', or 'uz'

Single Unit Rotate Test
*************************

This test takes a component and translates it along a specified axis. The function singleUnitRotateTest() takes the following arguments:

* index - The index of the component you would like to misalign, defaults to 0.
* min, max - The minimum and maximum rotational values you would like to test. An equal number of steps will be chosen between these values to be tested. Default to -1 deg and 1 deg, respectively.
   * min and max must be in units of length, see the section on :ref:`Astropy Units <units-top>`
* num - The number of different tests you would like to perform. The misalignment values will be equally spaced between the min and max values. Defaults to 10.
* axis - The dimension about which you would like to rotate. This argument should have a value of 1, 2, or 3 to rotate about the x, y, or z axis, respectively. Defaults to 1.
* plot - A boolean that tells the function if you would like to results to be automatically plotted. If True, a plot will be generated showing the full-width at half-max in the desired dimension as a function of misalignment value. Defaults to True.
* param - A string that tells the function the variable that you would like to plot the FWHM of. 
   * param has potential values: 'x', 'y', 'z', 'l', 'm', 'n', 'ux', 'uy', or 'uz'

Single Rotate Test
*************************

This test takes a component and translates it along a specified axis. The function singleRotateTest() takes the following arguments:

* index - The index of the component you would like to misalign, defaults to 0.
* min, max - The minimum and maximum rotational values you would like to test. An equal number of steps will be chosen between these values to be tested. Default to -1 deg and 1 deg, respectively.
   * min and max must be in units of length, see the section on :ref:`Astropy Units <units-top>`
* num - The number of different tests you would like to perform. The misalignment values will be equally spaced between the min and max values. Defaults to 10.
* ux,uy,uz - These argument describe the axis about which you would like to rotate. These arguments default to 1, 0, and 0, respectively.
* plot - A boolean that tells the function if you would like to results to be automatically plotted. If True, a plot will be generated showing the full-width at half-max in the desired dimension as a function of misalignment value. Defaults to True.
* param - A string that tells the function the variable that you would like to plot the FWHM of. 
   * param has potential values: 'x', 'y', 'z', 'l', 'm', 'n', 'ux', 'uy', or 'uz'

:ref:`Back to Top<instrument-top>`















