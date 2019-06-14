
Components
================================

Components are the essential building blocks of any simulation using PRTP.


Initializing a Component
------------------------------

Each Component in PRTP is governed by a class, where each class is located in a file. To import the necessary class, you need both the class name and the file name. For example, the Grating class is located within the Grating file. To import the class, you would use the syntax:

.. code-block:: python

   from prtp.Grating import Grating


Most standalone Components have identical class names and file names. However, Components which are part of a larger group are often grouped with these other Components in a single file. For example, the WolterOptic file stores every type of Wolter mirror. So if you wished to define just a Wolter Primary mirror, you could import the class using the syntax:

.. code-block:: python

   from prtp.WolterOptic import WolterPrimary


Once you've imported the necessary class, you just need to initialize an instance of the class. For a Grating, that could look like:

.. code-block:: python
   
   from prtp.Grating import Grating
   g = Grating(**kwargs)

In the above block, \**kwargs represents all of the necessary arguments. To see what arguments are necessary for a given Component type, see the page for that component.

Use an Existing Mission Component
----------------------------------

Included with PRTP are the components used in different missions that have already been created with all the proper parameters. All of these components are located within the prtp.Missions folder. To see all available missions and the Components each one includes, see the Mission page.

Each Mission Component is stored in a class within its mission's file. So to initialize a Mission Component, you must first initialize its class. Then, it is enough to call the class name with parenthesis to initialize the Component. For example, OGRE's Mirror Module can be accessed using the following lines of code:

.. code-block:: python

   from prtp.Missions.OGRE import OgreMirrorModule
   mirrormod = OgreMirrorModule()

Some Mission Components can accept additional arguments. For example, all Mission Sources can have the number of photons specified during initialization. To see the optional arguments for each Mission Component, see the Mission page.




