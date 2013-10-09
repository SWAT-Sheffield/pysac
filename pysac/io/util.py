import os
import numpy as np
from pysac.io import VACfile

def SAC_split_array(array,n0,n1,n2,axis_skip=0):
    """
    Split an array into the same order peices as SAC distribution.
    
    This is useful for splitting up a domain into bit to be distributed over mpi

    Notes
    -----
    This is only implemented for 3D    
    
    Parameters
    ----------
    array: np.ndarray
        The array to be split
    
    n0, n1, n2: int
        The number of splits along each axis
    
    axis_skip: int
        Skip the first n axes
    """
    zsplit = np.split(array,n2,axis=2+axis_skip)
    
    xsplit = []
    for z in zsplit:
        ysplit = []    
        ysplit.append(np.split(z,n1,axis=0+axis_skip))
        for y in ysplit:
            for x in y:
                xsplit += np.split(x,n0, axis=1+axis_skip)
    return np.array(xsplit)

def spilt_file(vac_data,n0,n1,n2,outfname):
    """
    Read in a VACdata class and save n0xn1xn2 VACfile classes following a 
    nameing template:
    
    Parameters
    ----------
    vac_data: io.VACdata
        The input file
    
    n0,n1,n2: int
        The number of axes splits
    
    outfname: str
        The output filename template (see below)
    
    Filename Template
    -----------------
    outfname should be a file in the form:
    /<path>/<filename>.<ext>
    
    The output will be:
    /<path>/<filename>_np<n0><n1><n2>_<00x>.<ext>
    """
    nx = vac_data.nx
    nx_out = [nx[0]/n0, nx[1]/n1, nx[2]/n2]
    
    header_out = vac_data.header
    header_out['nx'] = nx_out
    
    #Split arrays    
    x_split = SAC_split_array(vac_data.x,n0,n1,n2,axis_skip=1)
    w_split = SAC_split_array(vac_data.w,n0,n1,n2,axis_skip=1)
    
    #Filename processing:
    fileName, fileExtension = os.path.splitext(outfname)
    
    for n in range(n0*n1*n2):
        out = VACfile(fileName + '_np%02i%02i%02i_%03i'%(n0,n1,n2,n) + fileExtension, 
                         mode='w')
        
        out.header = header_out
        out.x = x_split[n]
        out.w = w_split[n]
        
        out.write_step()
        out.close()