.. _3D:

===========
3D Analysis
===========

Introduction
^^^^^^^^^^^^
The code in this submodule is used in *Mumford et. al. 2013* to create flux surfaces and perform the analysis.


Design and Code Layout
^^^^^^^^^^^^^^^^^^^^^^

The code uses MayaVi2 and the tvtk library which wrap vtk into a more 'pythonic' interface. There is some work being undertaken to port the code to pure vtk python bindings however the interface between numpy and vtk arrays is problematic and MayaVi provides a nice solution to this.

The analysis code is divided up into four major files and an animation script.