
Collimator Plate
===================

Collimator Plates are Flat Components, but they are specialized to handle complicated ways of eliminating photons which hit them. They can be used as pieces of a large Collimator, or individually to function as an Aperture.

Creating a Collimator Plate
-------------------------------

A Collimator Plate requires the following arguments:

* x,y,z - The spatial coordinates of the center of the plate. See :ref:`Flat Component <flat-component-definition>`
* nx,ny,nz - The components of the normal vector. See :ref:`Flat Component <flat-component-definition>`
* sx,sy,sz - The components of the surface vector. See :ref:`Flat Component <flat-component-definition>`
* l - The length of the Collimator Plate. This is the extent of the Component in the direction of the surface vector
   * If l is None, the length and width of the Collimator Plate will not be considered. That is, the Component will extent infinitely in both direction. If it is not None, it must be in units of length. See the section on Astropy Units.
* w - The width of the Collimator Plate. This is the extent of the Component in the direction of the cross product of the surface and normal vectors (sxn).
   * If w is not None, it must be in units of length. See the section on Astropy Units.
* collfunc - A function that defines how photons will be removed from the surface. More on this in a later section.

Moving a Collimator Plate
----------------------------

Collimator Plate objects inherit translate, rotate, and unitrotate from Flat Component, see the function usage :ref:`here <flat-component-motion>`

Collision Functions
--------------------

Collimator Plates have two ways to remove rays once they've been traced.
