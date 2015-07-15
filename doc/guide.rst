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
install it in 'develop' mode which enable you to edit the code in the directory
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
a python script or interpreter to open a SAC file you can now do::

    myfile = pysac.io.SACfile("mysacfile.out")

myfile now contains various routines to read the data in the file, as we shall 
see in more detail later.

Input/Output
^^^^^^^^^^^^

The first role of pySAC is to enable the reading and writing of SAC or VAC 
files in either the VAC style fortran unformatted binary file type or the 
newer HDF5 file structures.

The basic components of any VAC/SAC filetype are:

- x array: An array containing all coordinates at any point in the domain
- y array: An array containing all the computed variables
- header: A header, stored in some form consisting of at:
    - filehead - A not so descriptive file string
    - iteration number
    - physical time
    - number of dimensions
    - number of constant parameters in the equations
    - number of physcical varibles being solved for (in w array)
    - nx
    - the equation parameters
    - names for the varibles in the w array

pySAC implements a 'VACdata' class that will read either the VAC fortran type
files or HDF5 files containing VAC/SAC data. There is a specification of 
VACdata to SACdata for SAC files with magnetic field and pertubation and 
background varibles in the w array.

SACdata implemetes methods to extract primative variables such as temperature 
and pressure.

File Description
----------------

Header Contents
===============

The header of VAC / SAC files contains the same information irrespective of the
file type. When saving a file using pySAC the routines will automatically 
reformat the header to fit the file. The keys in the header and examples are 
listed below::

    header = {
                'filehead': '3Dmhd_33',
                'ndim': 3,
                'it': 0,
                't': 0,
                'nw': 13,
                'nx': [128,128,128],
                'varnames': ['z', 'x', 'y', 'h', 'm1', 'm2', 'm3', 'e', 'b1', 'b2', 'b3',
                             'eb', 'rhob', 'bg1', 'bg2', 'bg3', 'gamma', 'eta', 
                             'grav1', 'grav2', 'grav3'],
                'neqpar': 7,
                'eqpar': [1.66666667, 0.0, 1.0, -274.0, 0.0, 0.0, 0.0]
                }

Some extra parameters can be added for HDF5 files, namely: 'filedesc' which is 
a file desciption.

'varnames' should contain a label for each corrdinate dimension, a label for 
each varible in the w array, followed by each equation parameter. So ::

    len(header['varnames']) == nx + nw + neqpar

The above example is not consitent with this.

VAC File Description
====================
The unformatted FORTRAN binary files have the following format
for each time step stored a header of the following form preceeds the data:

- 'filehead' a string describing the file i.e. '2D_mhd22'
- 'params' a list containing:
    iteration number, physical time, ndim, number of equation params, number of vars in w array
    i.e. [10, 11.2541, 2, 7, 9], [int, float, int, int, int]
- nx coord dimensions e.g. [128,128,128]
- 'eqpar' - equation parameters, neqpars floats.
- varnames a list of strings nw long. Holds the names for all the w vars.
- x array
- w array

These are the defualt binary file type that VAC/SAC uses.

WARNING: These files are compiler and machine dependant, they are not portable
and should not be used over the far superior HDF5 files.
