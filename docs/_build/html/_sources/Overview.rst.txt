
Overview
================================

The Python Raytrace Package (PRTP) is a package of Python functions that simulate X-Ray instruments. It was designed for the McEntaffer Group at the Pennsylvania University. PRTP expands on the functionality of `PyXFocus <https://github.com/rallured/PyXFocus>`_, but with a more intuitive, modular design.

The foundation of PRTP lies in Components. These Components are Python objects which represent anything you would need in an x-ray instrument. They include mirrors, gratings, collimators, and more. Once they have been created with the proper position and orientation, tracing and reflecting photons can be done with just a few lines of code.

The whole package is optimized to allow straightforward instrument creation. Simply define your own components or pull components from an existing mission. From there, it is easy to test new designs or analyze alignment tolerances.

Installation
---------------------------

Currently, the only way to acquire PRTP is by downloading it from the github repository `here <https://github.com/bjmyers/prtp>`_.






