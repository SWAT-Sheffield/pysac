"""
A set of classes and routines for handling RAW SAC output in "unformatted FORTRAN"
files and also version 1 HDF5 files, before the adoption of the GDF layout
"""
import os
import struct

import h5py
import numpy as np


__all__ = ['VACfile', 'VAChdf5', 'VACdata', 'SACdata']

#==============================================================================
# File I/O Classes
#==============================================================================

class FortranFile(file):
    """
    File with methods for dealing with fortran unformatted data files

    Credit: Neil Martinsen-Burrell [via Enthought Mailing list]
    """

    def __init__(self,fname, mode='r', buf=0):
         """Open the file for writing, defaults to little endian."""
         file.__init__(self, fname, mode, buf)
         self.setEndian('<')

    def setEndian(self,c):
        """
        Set endian to big (c='>') or little (c='<') or native (c='@')

        Parameters
        ----------
        c : string
            The endian-ness to use when reading from this file.
        """
        if c == '<' or c == '>' or c =='@' or c == '=':
            self.ENDIAN = c
        else:
            raise ValueError('Cannot set endian-ness')

    def readString(self):
        """Read in a string with error checking"""
        l = struct.unpack(self.ENDIAN+'i',self.read(4))[0]
        str = self.read(l)
        if  struct.unpack(self.ENDIAN+'i',self.read(4))[0] != l:
            raise IOError('Error reading string from data file')
        return str

    def writeString(self,s):
        """
        Write a string

        Parameters
        ----------
        s : str
            The string to write
        """
        self.write(struct.pack(self.ENDIAN + 'i', len(s)))
        self.write(s.encode("latin-1"))
        self.write(struct.pack(self.ENDIAN + 'i', len(s)))

    def readReals(self, prec='d'):
        """Read in an array of reals (given precision) with error checking"""

        if prec not in ['d','f']:
             raise ValueError('Not an appropriate precision')

        l = struct.unpack(self.ENDIAN+'i',self.read(4))[0]
        data_str = self.read(l)
        len_real = struct.calcsize(prec)
        if l % len_real != 0:
            raise IOError('Error reading array of reals from data file')
        num = l/len_real
        reals = struct.unpack(self.ENDIAN+str(num)+prec,data_str)
        if struct.unpack(self.ENDIAN+'i',self.read(4))[0] != l:
            raise IOError('Error reading array of reals from data file')
        return list(reals)

    def writeReals(self, reals, prec='d'):
        """
        Write an array of floats in given precision

        Parameters
        ----------
        reals : array
            Data to write
        prec : str
            Character code for the precision to use in writing
        """
        reals = np.asarray(reals)
        if prec not in ['d','f']:
            raise ValueError('Not an appropriate precision')

        self.write(struct.pack(self.ENDIAN + 'i', reals.nbytes))
        self.write(struct.pack(self.ENDIAN + prec * reals.size,
                               *reals.flatten(order='F')))
        self.write(struct.pack(self.ENDIAN + 'i', reals.nbytes))

    def readInts(self):
        """Read in an array of integers with error checking"""
        l = struct.unpack('i',self.read(4))[0]
        data_str = self.read(l)
        len_int = struct.calcsize('i')
        if l % len_int != 0:
            raise IOError('Error reading array of integers from data file')
        num = l/len_int
        ints = struct.unpack(str(num)+'i',data_str)
        if struct.unpack(self.ENDIAN+'i',self.read(4))[0] != l:
            raise IOError('Error reading array of integers from data file')
        return list(ints)

    def writeInts(self, ints):
        """
        Write an array of integers in given precision

        Parameters
        ----------
        ints : array
            Data to write
        """
        self.write(struct.pack(self.ENDIAN + 'i', struct.calcsize('i') * len(ints)))
        self.write(struct.pack(self.ENDIAN + 'i' * len(ints), *ints))
        self.write(struct.pack(self.ENDIAN + 'i', struct.calcsize('i') * len(ints)))

    def readRecord(self):
         """Read a single fortran record"""
         l = struct.unpack(self.ENDIAN+'i',self.read(4))[0]
         data_str = self.read(l)
         # check length
         if len(data_str) != l:
             raise IOError('Didn''t read enough data')
         check = self.read(4)
         if len(check) != 4:
             raise IOError('Didn''t read enough data')
         if struct.unpack(self.ENDIAN+'i',check)[0] != l:
             raise IOError('Error reading record from data file')
         return data_str

    def readParams(self,prec='d'):
        """Reads the Params line which is a mix of Ints and Reals"""
        #Check that prec is spec'd proper
        if prec not in ['d','f']:
            raise ValueError('Not an appropriate precision')
        #read in line
        data_str = self.readRecord()
        pars = struct.unpack(self.ENDIAN+'idiii', data_str)
        return list(pars)

    def writeParams(self, params):
        self.write(struct.pack(self.ENDIAN + 'i', 24))
        self.write(struct.pack(self.ENDIAN + 'idiii', *params))
        self.write(struct.pack(self.ENDIAN + 'i', 24))

class VACfile():
    def __init__(self, fname, mode='r', buf=0):
        """
        Base input class for VAC Unformatted binary files.
        Based on FortranFile has been modified to read VAC / SAC output files.

        Parameters
        ----------
        fname:    string
            Input Filename

        mode: {'r' | 'w'}
            I/O mode (only 'r' is fully supported)

        buf: int
            underlying I/O buffer size
        Returns
        -------
        Reads a iteration into the following structure:
        
        file.header: -Dictionary containging
                     -filehead: string at begging of file
                     -params: Iteration Parameters, it, t, ndim, neqpar, nw
                     -nx: Size of cordinate array [list]
                     -eqpar: eqpar\_ parameters [list]
                     -varnames: list containg varible names for dimensions, nw and eqpar?
        file.w : w array from file which is [params,[nx]] in size
        file.w\_: dict containing the {varname:index} pairs for the w array
        file.x : x array from file which is [ndim,[nx]] in size
        
        """
        #Do FORTRAN read init, set Endian for VAC/SAC files
        self.file = FortranFile(fname,mode,buf)
        self.file.setEndian('<')

        if mode == 'r':
            self.process_step()
            self.recordsize = self.file.tell()
            self.num_records = os.stat(fname).st_size / self.recordsize

            #Find out first and last time values
            self.t_start = self.header['params'][1]
            self.read_timestep(self.num_records)
            self.t_end = self.header['params'][1]
            self.header['final t'] = self.t_end
            self.read_timestep(1)

            print "File is %i Records Long"%self.num_records

        elif mode == 'w':
            #Set up varibles to be filled and then written
            self.header = []
        else:
            raise ValueError("mode must be 'r' or 'w'")

    def read_timestep(self,i):
        self.file.seek(int(i-1) * self.recordsize)
        self.process_step()

    def process_step(self):
        """
        Does the raw file processing for each timestep

        Sets up the header and reads and reshapes the arrays
        """
        self.header = {}
        self.header['filehead'] = self.file.readString()
        self.header['params'] = self.file.readParams()
        #params is: it, t, ndim, neqpar, nw
        self.header['it'] = self.header['params'][0]
        self.header['t'] = self.header['params'][1]
        self.header['ndim'] = self.header['params'][2]
        self.header['neqpar'] = self.header['params'][3]
        self.header['nw'] = self.header['params'][4]
        self.header['nx'] = self.file.readInts()
        self.header['eqpar'] = self.file.readReals()
        self.header['varnames'] = self.file.readRecord().split()
        
        # Reading in the FORTRAN arrays into Python can be achieved in a
        # number of different ways.
        # The arrays on disk are in the following order:
        # [z, x, y, n]
        # Where n is the coordinate in the x array and the varible in the case
        # of the w array.
        # Note: the w array is actually saved as nw * [z, x, y] arrays rather 
        # than one 4D array, unlike the x array.
        
        # X ARRAY
        # Reading in the x array in the same order that it is in on disk:
        self.x = self.file.readReals()
        s = self.header['nx'] + [self.header['ndim']]
        self.x = np.reshape(self.x, s, order='F')
        # This returns [z, x, y, 3] for a 3D array, where the indicies of the 
        # last dimension are 0 - z, 1 - x, 2 - y.
        
# This code might read x in with the index dimension first, it is also possible
# to achieve that using reshape. But it is safer to keep it in disk order.
#        self.x = np.zeros([self.header['ndim']] + self.header['nx'])
#        for i in range(0,self.header['ndim']):
#            self.x[i] = np.reshape(self.file.readReals(), self.header['nx'], order='F')
        

        # The w array is read in as a sequence of 3D arrays, so the index 
        # dimension is put first.
        self.w = np.zeros([self.header['params'][-1]]+self.header['nx'], order='F')
        for i in xrange(0, self.header['params'][-1]):
            self.w[i] = np.reshape(self.file.readReals(), self.header['nx'], order='F')
        self.w_ = {}
        ndim = self.header['params'][2]
        nw = self.header['params'][-1]
        #find h in varnames (to prevent x y h bug in 3D file)
        index = next((i for i in xrange(len(self.header['varnames']))
                    if not(self.header['varnames'][i] in ["x","y","z"])),ndim)
        for i,name in enumerate(self.header['varnames'][index:nw+index]):
            self.w_.update({name:i})

    def _write_header(self):
        """ Write header information """

        #Validate Header
        if not all([t in self.header for t in ['filehead', 'nx', 'eqpar']]):
            raise ValueError('Invalid header for writing')

        if not ('params' in self.header or
                all([t in self.header for t in
                                        ['it', 't', 'ndim', 'neqpar', 'nw']])):
            raise ValueError('Invalid header for writing')

        #pad the header to 79 characters as the distribute script needs it like
        #that even though SAC dosen't
        self.file.writeString(self.header['filehead'])#.ljust(79))

        if 'params' in self.header:
            params = self.header['params']

        else:
            params = [self.header['it'], self.header['t'], self.header['ndim'],
                      self.header['neqpar'], self.header['nw']]

        self.file.writeParams(params)
        self.file.writeInts(self.header['nx'])
        self.file.writeReals(self.header['eqpar'])
        self.file.writeString(" ".join(self.header['varnames']))

    def _write_data(self):
        """ Save arrays into unformatted fortran file """
        self.file.writeReals(self.x, prec='d')
        for w in self.w:
            self.file.writeReals(w, prec='d')

    def write_step(self):
        #Make sure you are saving the correct size data
        assert tuple(self.header['nx']) == self.w[0].shape
        assert tuple(self.header['nx']) == self.x[0].shape

        self._write_header()
        self._write_data()

    def close(self):
        self.file.close()

class VAChdf5():
    def __init__(self, filename, mode='r'):
        """
        Based on FortranFile has been modified to read VAC / SAC HDF5 files.

        Reads a iteration into the following structure:
        
        file.header: 
        -Dictionary containging
        -filehead: string at begging of file
        -params: Iteration Parameters, it, t, ndim, neqpar, nw
        -nx: Size of cordinate array [list]
        -eqpar: eqpar\_ parameters [list]
        -varnames: list containg varible names for dimensions, nw and eqpar?
        file.w : w array from file which is [params,[nx]] in size
        file.w\_: dict containing the {varname:index} pairs for the w array
        file.x : x array from file which is [ndim,[nx]] in size

        Also creates HDF5 specific attributes:
        
        file.sac_group - Holds the x and time_group attributes.
        file.time_group - Holds the series of w arrays.

        Largely the HDF5 file is designed so the functionality mimics the VAC
        binary file, i.e. all the vars are still in the W array etc.
        """
        self.mode = mode
        self.h5file = h5py.File(filename,mode)

        if mode == 'r':
            #Open top level group
            if not("SACdata" in self.h5file.keys()):
                print """Are you sure this is a proper SAC HDF5 file?
                    Opening first group."""
                self.sac_group = self.h5file[self.h5file.keys()[0]]
            else:
                self.sac_group = self.h5file["SACdata"]

            self.x = self.sac_group['x']

            self.header = dict(self.sac_group.attrs)
            self.header['neqpar'] = int(self.header['neqpar'])
            self.header['ndim'] = int(self.header['ndim'])
            try:
                self.header['filehead'] = self.h5file.attrs['filehead'][0]
            except:
                pass

            self.time_group = self.sac_group['wseries']
            self.header.update(dict(self.time_group.attrs))
            self.header['varnames'] = self.header['varnames'][0].split()
            self.read_timestep(0)
            self.t_start = self.header['t']
            if 'final t' in self.sac_group.attrs:
                self.header['final t'] = self.sac_group.attrs['final t']
                self.t_end = self.sac_group.attrs['final t']
            else:
                self.t_end = self.time_group.items()[-1][1].attrs['t']#[0] ##don't know
                self.header['final t'] = self.t_end

            self.num_records = len(self.time_group.items())

        elif mode =='w':
            self._has_header = False
            #Create sacdata group
            self.sac_group = self.h5file.create_group("SACdata")
            #create wseries group
            self.time_group = self.sac_group.create_group("wseries")
        else:
            raise ValueError("mode must be 'r' or 'w'")

    def read_timestep(self,i):
        wstepname = self.time_group.keys()[i]
        self.header['it'] = int(self.time_group[wstepname].attrs['it'])
        self.header['t'] = float(self.time_group[wstepname].attrs['t'])
        #to maintain backwards compatibility with VACfile
        self.header['params'] = [self.header['it'], self.header['t'],
                                self.header['ndim'], self.header['neqpar'],
                                self.header['nw']]

        self.w = self.time_group[wstepname]
        self.w_ = {}
        self.w_dict = {}
        index = next((i for i in xrange(len(self.header['varnames'])) if not(
            self.header['varnames'][i] in ["x","y","z"])),self.header['ndim'])

        for i,name in enumerate(
                    self.header['varnames'][index:self.header['nw'] + index]):
            self.w_.update({name:i})
            self.w_dict.update({name:self.w[i]})

    def _write_header(self):
        """ Populate top level attributes with the file meta data """
        if self.mode == 'r':
            raise Exception("This is only vaild if mode = 'r'")

        #Write file level attributes
        self.h5file.attrs.create('filehead', self.header['filehead'])
        if 'filedesc' in self.header:
            desc = self.header['filedesc']
        else:
            desc = 'This is a SAC HDF5 file written by pySAC'
        self.h5file.attrs.create('filedesc', desc)

        #populate the SAC group
        self.sac_group.attrs.create('eqpar', self.header['eqpar'])
        self.sac_group.attrs.create('ndim', self.header['ndim'])
        self.sac_group.attrs.create('neqpar', self.header['neqpar'])
        self.sac_group.attrs.create('nx', self.header['nx'])

        #Populate the wseries group
        self.time_group.attrs.create('nw', self.header['nw'])
        #This is saved in a list to match FORTRAN behavior
        self.time_group.attrs.create('varnames', [" ".join(self.header['varnames'])])

        #write x array
        self.sac_group.create_dataset('x', data=self.x)

        self._has_header = True

    def write_step(self):
        """ Save step data into hdf5 file """
        if self.mode == 'r':
            raise Exception("This is only vaild if mode = 'r'")

        if not self._has_header:
            self._write_header()

        dset = self.time_group.create_dataset('w%05i'%self.header['it'], data=self.w)

        dset.attrs.create('it', self.header['it'])
        dset.attrs.create('t', self.header['t'])

    def close(self):
        if self.mode == 'w':
            #Write out final info to sacgroup
            self.sac_group.attrs.create('final t', self.header['t'])
            self.sac_group.attrs.create('nt', self.header['it'])
        self.h5file.close()


#==============================================================================
# VAC / SAC data classes for UI and Output (hopefully)
#==============================================================================

class VACdata(object):
    """ This is a file type independant class that should expose a VACfile or
    VACHDF5 file so it is transparent to the end user. """

    def _get_file_type(self, filename, filetype='auto'):
        """ This method determines the filetype based on user input or the
        file extension"""

        filePrefix, fileExt = os.path.splitext(filename)
        if fileExt == '.h5':
            ofiletype = 'hdf5'
        elif fileExt in ['.ini', '.out']:
            ofiletype ='fort'
        else:
            if filetype == 'auto':
                raise TypeError(
                        "File type can not be automatically determined")
            else:
                if filetype in ['hdf5', 'fort']:
                    ofiletype = filetype
                else:
                    raise ValueError(
"Specified filetype is not valid. Filetype should be one of { 'hdf5' | 'fort' }")
        return ofiletype

    def __init__(self, filename, filetype='auto'):
        """
        Create a VACdata class.

        Parameters
        ----------
        filename: str

        filetype: str {'auto' | 'fort' | 'hdf5' }
        """

        #Detect filetype and open file for reading.
        filetype = self._get_file_type(filename, filetype)

        if filetype == 'hdf5':
            self.file = VAChdf5(filename)

        elif filetype == 'fort':
            self.file = VACfile(filename, mode='r')

    @property
    def w(self):
        return self.file.w

    @property
    def x(self):
        return self.file.x

    @property
    def w_(self):
        return self.file.w_

    @property
    def w_dict(self):
        return self.file.w_dict

    @property
    def header(self):
        return self.file.header

    @property
    def t_start(self):
        return self.file.t_start

    @property
    def t_end(self):
        return self.file.t_end

    @property
    def num_records(self):
       return self.file.num_records

    def read_timestep(self, i):
        """
        Read in the specified time step

        Parameters
        ----------
        i: int
            Time Step number
        """
        self.file.read_timestep(i)

class SACdata(VACdata):
    """
    This adds specifications to VACdata designed for SAC simulations in 2D or
    3D with magnetic field.

    This adds the background and pertubation varibles into a new w_sac dict.
    """

    def __init__(self, filename, filetype='auto'):
        VACdata.__init__(self, filename, filetype=filetype)
        self.update_w_sac()

    def update_w_sac(self):
        """
        This method creates the w_sac dictionary for the current timestep.
        """
        self.w_sac = {}
        if self.header['ndim'] == 2:
            self.w_sac.update({'rho':self.w[self.w_["h"]] +
                                                    self.w[self.w_["rhob"]]})
            self.w_sac.update({'v1':self.w[self.w_["m1"]] / self.w_sac['rho']})
            self.w_sac.update({'v2':self.w[self.w_["m2"]] / self.w_sac['rho']})
            self.w_sac.update({'e':self.w[self.w_["e"]] +
                                                        self.w[self.w_["eb"]]})
            self.w_sac.update({'b1':self.w[self.w_["b1"]] +
                                                    self.w[self.w_["bg1"]]})
            self.w_sac.update({'b2':self.w[self.w_["b2"]] +
                                                    self.w[self.w_["bg2"]]})

        if self.header['ndim'] == 3:
            self.w_sac.update({'rho':self.w[self.w_["h"]] +
                                                    self.w[self.w_["rhob"]]})
            self.w_sac.update({'v1':self.w[self.w_["m1"]] / self.w_sac['rho']})
            self.w_sac.update({'v2':self.w[self.w_["m2"]] / self.w_sac['rho']})
            self.w_sac.update({'v3':self.w[self.w_["m3"]] / self.w_sac['rho']})
            self.w_sac.update({'e':self.w[self.w_["e"]] +
                                                        self.w[self.w_["eb"]]})
            self.w_sac.update({'b1':self.w[self.w_["b1"]] +
                                                    self.w[self.w_["bg1"]]})
            self.w_sac.update({'b2':self.w[self.w_["b2"]] +
                                                    self.w[self.w_["bg2"]]})
            self.w_sac.update({'b3':self.w[self.w_["b3"]] +
                                                    self.w[self.w_["bg3"]]})

    def get_w_yt(self):
        self.w_yt = {}
        if self.header['ndim'] == 3:
            self.w_yt.update({'Bx':self.w_sac['b2'],
                              'By':self.w_sac['b3'],
                              'Bz':self.w_sac['b1']})

            self.w_yt.update({'x-velocity':self.w_sac['v2'],
                              'y-velocity':self.w_sac['v3'],
                              'z-velocity':self.w_sac['v1']})

            self.w_yt.update({'Density':self.w_sac['rho'],
                              'e':self.w_sac['e']})

        if self.header['ndim'] == 2:
            raise ValueError("Doesn't support 2D")

        return self.w_yt

    def read_timestep(self,i):
        VACdata.read_timestep(self,i)
        self.update_w_sac()

    def convert_B(self):
        """
        This function corrects for the scaling of the magentic field units.

        It will convert the magnetic field into Tesla for the current time
        step.

        WARNING: The conservative variable calculations are in SAC scaled
        magnetic field units, this conversion should be run after accessing any
        calculatuions involving the magnetic field
        """
        mu = 1.25663706e-6
        if self.header['ndim'] == 2:
            self.w_sac['b1'] *= np.sqrt(mu)
            self.w_sac['b2'] *= np.sqrt(mu)
        if self.header['ndim'] == 3:
            self.w_sac['b1'] *= np.sqrt(mu)
            self.w_sac['b2'] *= np.sqrt(mu)
            self.w_sac['b3'] *= np.sqrt(mu)

    def get_thermalp(self,beta=False):
        """Calculate Thermal pressure from varibles """
        if self.header['ndim'] == 3:
            #raise NotImplementedError("This Dosen't work for 3D yet, go fix")
            g1 = (self.header['eqpar'][0]-1)
            kp = (self.w_sac['rho'] * (self.w_sac['v1']**2 + self.w_sac['v2']**2 + self.w_sac['v3']**2))/2.
            mp = (self.w_sac['b1']**2 + self.w_sac['b2']**2 + self.w_sac['b3']**2) / 2.
            p = g1 * (self.w_sac['e'] - kp - mp)
            #p = (\gamma -1) ( e - \rho v^2/2 - B^2/2)
        else:
            g1 = (self.header['eqpar'][0]-1)
            kp = (self.w_sac['rho'] * (self.w_sac['v1']**2 + self.w_sac['v2']**2))/2.
            mp = (self.w_sac['b1']**2 + self.w_sac['b2']**2) / 2.
            p = g1 * (self.w_sac['e'] - kp - mp)

        if beta:
            return p, mp
        else:
            return p

    def get_bgp(self):
        print "WARNING: Background Pressure will not work if inital conditions are not V=0"
        if self.header['ndim'] == 3:
            #raise NotImplementedError("This Dosen't work for 3D yet, go fix")
            g1 = (self.header['eqpar'][0]-1)
            kp = 0.0#(self.w[self.w_["rhob"]] * (self.w_sac['v1']**2 + self.w_sac['v2']**2 + self.w_sac['v3']**2))/2.
            mp = (self.w[self.w_["bg1"]]**2 + self.w[self.w_["bg2"]]**2 + self.w[self.w_["bg3"]]**2) / 2.
            p = g1 * (self.w[self.w_["eb"]] - kp - mp)
            #p = (\gamma -1) ( e - \rho v^2/2 - B^2/2)
        else:
            g1 = (self.header['eqpar'][0]-1)
            kp = 0.0#(self.w[self.w_["rhob"]] * (self.w_sac['v1']**2 + self.w_sac['v2']**2))/2.
            mp = (self.w[self.w_["bg1"]]**2 + self.w[self.w_["bg2"]]**2) / 2.
            p = g1 * (self.w[self.w_["eb"]] - kp - mp)
        return p

    def get_total_p(self):
        if self.header['ndim'] == 3:
           gamma = self.header['eqpar'][0]

           vtot2 = (self.w_sac['v1']**2 + self.w_sac['v2']**2 + self.w_sac['v3']**2)
           therm = self.w[self.w_["e"]] - (self.w_sac["rho"] * vtot2) / 2.

           Bpert = self.w[self.w_['b1']] + self.w[self.w_['b2']] + self.w[self.w_['b3']]
           Bpert2 = self.w[self.w_['b1']]**2 + self.w[self.w_['b2']]**2 + self.w[self.w_['b3']]**2
           Bback = self.w[self.w_['bg1']] + self.w[self.w_['bg2']] + self.w[self.w_['bg3']]
           mag = Bback * Bpert + (Bpert2 / 2.)

           return (gamma - 1) * therm - (gamma - 2) * mag
        else:
            raise NotImplementedError("This Dosen't work for 2D yet, go fix")

    def get_temp(self,p=None):
        if not(p):
            p = self.get_thermalp()
        T = (p * 1.2) / (8.3e3 * self.w_sac['rho'])
        return T

    def get_bgtemp(self):
        print "WARNING: Background Temprature will not work if inital conditions are not V=0"
        if self.header['ndim'] == 3:
            kp = 0.0#(self.w[self.w_["rhob"]] * (self.w_sac['v1']**2 + self.w_sac['v2']**2 + self.w_sac['v3']**2))/2.
            mp = (self.w[self.w_["bg1"]]**2 + self.w[self.w_["bg2"]]**2 + self.w[self.w_["bg3"]]**2) / 2.
            T = self.w[self.w_["eb"]] - kp - mp
        else:
            kp = 0.0#(self.w[self.w_["rhob"]] * (self.w_sac['v1']**2 + self.w_sac['v2']**2))/2.
            mp = (self.w[self.w_["bg1"]]**2 + self.w[self.w_["bg2"]]**2) / 2.
            T = self.w[self.w_["eb"]] - kp - mp
        return T

    def get_va(self):
        return (np.sqrt(self.w_sac['b1']**2 + self.w_sac['b2']**2
                        + self.w_sac['b3']**2) / np.sqrt(self.w_sac['rho']))
        #return (abs(self.w_sac['b1']) + abs(self.w_sac['b2']) + abs(self.w_sac['b3'])) / sqrt(self.w_sac['rho'])

    def get_cs(self,p=None):
        if not p:
            p = self.get_thermalp()
        g1 = self.header['eqpar'][0]
        return np.sqrt((g1 * p) / self.w_sac['rho'])
