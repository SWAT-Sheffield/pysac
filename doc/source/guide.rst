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

pySAC implements a `VACdata' class that will read either the VAC fortran type
files or HDF5 files containing VAC/SAC data. There is a specification of 
VACdata to SACdata for SAC files with magnetic field and pertubation and 
background varibles in the w array.

SACdata implemetes methods to extract primative variables such as temperature 
and pressure.

VAC File Description
--------------------
The unformatted FORTRAN binary files have the following format
for each time step stored a header of the following form preceeds the data:

- 'filehead' a string describing the file i.e. '2D_mhd22'
- 'params' a list containing:
    iteration number, physical time, ndim, number of equation params, number of vars in w array
    i.e. [10, 11.2541, 2, 7, 9], [int, float, int, int, int]
- nx number of dimensions (int)
- 'eqpar' - equation parameters, neqpars floats.
- varnames a list of strings nw long. Holds the names for all the w vars.
- x array
- w array

Input
=====

In this section we shall look at reading files using pySAC

VAC FORTRAN Files
-----------------

These are the defualt binary file type that VAC/SAC uses.

WARNING: These files are compiler and machine dependant, they are not portable
and should not be used over the far superior HDF5 files.

SAC HDF5 Files
--------------

This file type has been added to SAC to make the output standard and portable,
also to enable parallel I/O.

Output
======

Output in pySAC.io is done by writing out the current state of the VACdata 
object. To create a new file with data from elsewhere you would create a VACdata
object with mode='w' and then assign data and the header:
    
    myfile = sacio.VACdata("myoutfile.h5")
    myfile.header = header
    myfile.w = w_arr
    myfile.x = x_arr

Each time step can then be written by a call to write_step:
    
    myfile.write_step()

remember to close the file when you are done with:
    
    myfile.close()

The output routines will automatically determine the file type.

It is also possible to save out to a different file, or file_type by first 
reading in a file:
    
    myfile = sacio.VACdata("myinfile.h5")

then calling init_file() and write_step() for each iteration in the file:
    
    myfile.init_file("myoutfile.out")
    for i in range(num_records):
        myfile.read_timestep(i)
        myfile.write_step()
    myfile.close()

