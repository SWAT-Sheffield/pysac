===========
pySAC Guide
===========

pySAC is a python package that is designed to enable analysis and plotting of
VAC (Verstatile Advection Code) and SAC (Sheffield Advanced Code) as well as 
smaug, the GPGPU port of SAC.

Installation
^^^^^^^^^^^^

pySAC is contained inside a DVS version control system called git, with a 
hosted repository on BitBucket_. To obtain the latest version of the source
code you can clone the git repository thus::

    git clone git@bitbucket.org:swatsheffield/pysac.git

which you can then install using the included python setup script. There are
two options at this point, you can install the code permantly or you can 
install it in `develop' mode which enable you to edit the code in the directory
in which you have just cloned it.
To install in devlop mode do the following::
    
    python setup.py develop

which, under linux you will need to run as root.

.. _BitBucket: https://bitbucket.org/swatsheffield/pysac/

The following python packages are required to use pySAC:

- matplotlib >= 1.1
- h5py
- numpy
- scipy


Basic Usage
-----------

Once installed pySAC can be used by importing the component you wish to use in
a python script or interpreter::

    >>> import pysac.io

to open a SAC file you can no do::

    >>> myfile = pysac.io.SACdata("mysacfile.h5")

myfile now contains various routines to read the data in the file, as we shall 
see in more detail later.

Input/Output
------------

The first role of pySAC is to enable the reading and writing of SAC or VAC 
files in either the VAC style fortran unformatted binary file type or the 
newer HDF5 file structures.

The basic components of any VAC/SAC filetype are:

- x array: An array containing all coordinates at any point in the domain
- y array: An array containing all the computed variables
- header: A header, stored in some form consisting of at:
    - iteration number
    - physical time
    - number of dimensions
    - number of constant parameters in the equations
    - number of physcical varibles being solved for (in w array)
    - nx
    - the equation parameters
    - names for the varibles in the w array

pySAC implements a `VACdata' class that will read either the VAC fortran type
files or HDF5 files containing VAC/SAC data. There is a specification of 
VACdata to SACdata for SAC files with magnetic field and pertubation and 
background varibles in the w array.

SACdata implemetes methods to extract primative variables such as temperature 
and pressure.



