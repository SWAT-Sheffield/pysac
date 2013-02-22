# -*- coding: utf-8 -*-
"""
:Created on: Thu May 17 16:56:35 2012

:author: Stuart Mumford

Library for doing 2D slices and animations of SAC / VAC output
"""
from __future__ import absolute_import, division

import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from mpl_toolkits.axes_grid1 import make_axes_locatable
#from scipy import ndimage
from pysac.plot.CustomColourmaps import *

class FixedCentre(mpl.colors.Normalize):
    """
    Normalise with a Fixed Centre
    """
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale(result)
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result.fill(0.5)   # Or should it be all masked?  Or 0.5?
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)
            # ma division is very slow; we can take a shortcut
            resdat = result.data
            resdat[(result.data == 0.).nonzero()] = 0.0
            resdat[(result.data < 0.).nonzero()] /= (-vmin * 2.)
            resdat[(result.data > 0.).nonzero()] /= (vmax * 2.)
            result = np.ma.array(resdat+0.5, mask=result.mask, copy=False)
        if is_scalar:
            result = result[0]
        return result



class BaseArtist():
    """ Class to hold drawing routines for any type of image, line etc. """
    def plot():
        raise NotImplementedError("This needs to be subclassed")
    
    def update():
        raise NotImplementedError("This needs to be subclassed")
        

class Image(BaseArtist):
    """ Plots an imageartist on subplot """
    pass
            
class Subplot():
#        #create the subplot properties dictionary
#        subaxprops = {'subkey':self.sub_key, 'title':title,
#                      'fieldlines':fieldlines, 'colorbar':colorbar,
#                      'tbar':tbar, 'labels':labels,'nonuniform':nonuniform,
#                      'arr_slice':arr_slice,'cmap':cmap,'norm':norm}   
    def __init__(self,subkey,datakey):
        self.subkey = subkey
        self.datakey = datakey
        self.artists = []
    def set_cmap(self,cmap):
        if (cmap and not isinstance(cmap,mpl.colors.Colormap)):
            #Then it should be a string
            if cmap == 'jetblack':
                self.cmap = jetblack
            elif cmap in mpl.cm.cmap_d.keys():
                self.cmap = mpl.cm.get_cmap(name = cmap)
            else:
                raise Exception("""This is not a valid colormap, specify:
                    None for default behaviour or
                    an instance of mpl.colors.Colormap
                    or a string which is a default mpl.cm
                    or 'jetblack' """)
        else:
            self.cmap = cmap
    
    def set_norm(self,norm):
        if (norm and not isinstance(norm,mpl.colors.Normalize)):
            if norm == 'fixedcentre':
                self.norm = FixedCentre()
            else:
                raise Exception("""This is not a valid mpl.colors.Normalize
                instance, please specifiy an instance of mpl.colors.Normalize
                or 'fixedcentre for SACplot.FixedCentre()""")
        else:
            self.norm = norm

        
    

#Is this class going to have to be a subclasss of mpl.animation.timedanim???
class SACplot():
    def __init__(self, sacfile,sup_title='',figsize=(12,4)):
        """A Visualistion interface for SAC data files and simulations.
        
        This interface will when provided with a SACfile instance can be used 
        to specify and create either a plot or an animation of the data.
        
        All primitive or conservative varibles created in the SAC simulation
        """
        self.f = sacfile #TODO: Make this a filename??!
        self.fig = plt.figure(figsize=figsize)
        self.fig.suptitle(sup_title,fontsize=16)
        #Create a list of subplot objects
        self.subplots = []
        self.sub_key = 1        
        
        #Fix so only generate fieldlines once per frame
        self.fieldseeds = False
        self.line_colour = 'orange'
        #This prevents multiple file reads
        self.currentrecord = None
        self.explots = []

    def grid_spec(self,nc,nr,**kwargs):
        """ Set Row / Coloumn layout of subplots and subplotparams """
        self.nc = nc
        self.nr = nr
        
        if kwargs:
            self.fig.subplots_adjust(kwargs)
        else:
            self.fig.subplots_adjust(top=0.88,wspace=0.2,left=0.05,
                                      right=0.93,hspace=0.4,bottom=0.05)
    def add_tbar(self):
        self.HAS_tbar = True
        bax = self.fig.add_axes([0.4,0.93,0.2,0.01])
        self.tbar = bax.barh(0,self.f.header['t'],color='k') #TODO: Anim FIX!
        bax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=100.0))
        bax.yaxis.set_major_locator(mpl.ticker.NullLocator())
        bax.set_xlim((self.f.t_start,self.f.t_end))
        bax.set_ylabel(r"$t=$",rotation='horizontal')


    def add_img(self,datakey,title='',fieldlines=False,colorbar=False,labels=[],
                arr_slice=np.s_[...],cmap=None,norm=None):
    
        #TODO: This can be removed in favor of auto subplot counting each time
        # this routine is called.
        #HOWEVER: Need a way to determine number of coloumns?!
        if self.sub_key > self.nc * self.nr:
            raise Exception("grid_spec wrong or too many plots added!")
            
        #Create a subplot object
        subplot = Subplot(self.sub_key,datakey)
        subplot.colorbar = colorbar
        subplot.labels = list(labels)
        subplot.arr_slice = arr_slice
        subplot.title = title
        

        subplot.set_cmap(cmap)
        
        subplot.set_norm(norm)
        
        if isinstance(fieldlines,list):
            #if self.fieldseeds is already defined just say NO
            if isinstance(self.fieldseeds,list):
                raise ValueError("You can only specify fieldline parameters in the first frame")
            #Where fieldlines is a list of parameters make sure its 3 pars long
            if len(fieldlines) != 3:
                raise Exception("""fieldlines should be list: [ax,denisty,dS] 
                                or boolean""" )
            ax = fieldlines[0]
            density = fieldlines[1]
            dS = fieldlines[2]
            #make the lines
            self.fieldseeds = self.seed_generator(ax,density)
            self.dS = dS
                    
        elif isinstance(fieldlines,bool):
           pass #Store in subaxprobs and be merry
           
        elif not(isinstance(fieldlines,bool)) or not(isinstance(fieldlines,list)):
            raise TypeError("fieldlines argument must be Boolean or List")
            
        else:
            #Catch all
            raise Exception("Something is wrong with fieldlines argument")
        subplot.fieldlines = fieldlines

        self.subplots.append(subplot)
        self.sub_key += 1
    
    def _read_data(self,keys,i):
        """Read and process a step of data."""
        if self.currentrecord != i:
            #Read in requested step
            self.f.readrecord(i)
        self.currentrecord =  i
        if type(keys) is list:
            all_data = []
            for key in keys:
                all_data.append(self._read_key(key))
            return all_data
        else:
            return self._read_key(keys)
            
    def _read_key(self,key):
        #If it's in w_sac then easy:
        if self.f.w_sac.has_key(key):
            data = self.f.w_sac[key]
        elif self.f.w_.has_key(key):
            data = self.f.w[self.f.w_[key]]
        #else do it yourself:
        elif key == 'pressure':
            data = self.f.get_thermalp()
        elif key == 'temp':
            data = self.f.get_temp()
        elif key == 'delta_p':
            bgp = self.f.get_bgp()
            data = (self.f.get_thermalp() - bgp)/bgp
        elif key == 'delta_t':
            bgT = self.f.get_bgtemp()
            data = (self.f.get_temp() - bgT)/bgT
        elif key == 'v_total':
            data = np.sqrt(self.f.w_sac['v2']**2 + self.f.w_sac['v3']**2)
        elif key == 'b_total':
            data = np.sqrt(self.f.w_sac['b1']**2 + self.f.w_sac['b2']**2)
        elif key == 'delta_rho':
            data = self.f.w[self.f.w_['h']] - self.f.w[self.f.w_['rhob']]
        elif key == 'va':
            data = self.f.get_va()
            data = np.abs(data)
        elif key == 'cs':
            data = self.f.get_cs()
        else:
            #Seriously you should have raised two errors before this!
            raise Exception("Please make sure your keys are right!")
        return data
        
    def plot(self,i):
        """Reads data, and plots requested subplots. Last routine to call."""
        
        for subplot in self.subplots:
            self.get_subplot_extent(subplot)
            subplot.data = self._read_data(subplot.datakey,i)
            subplot.data = subplot.data[subplot.arr_slice]
            subplot.t = [self.f.header['t']]
#            import pdb; pdb.set_trace()
            if self.fieldseeds and subplot.fieldlines:
                #Make and store fieldlines for this datapoint
                ##This almost certainly could be done better::
                if self.f.ndim == 2:
                    B2 = self._read_data('b1',i)
                    B1 = self._read_data('b2',i)
                else:
                    B2 = self._read_data('b%i'%(subplot.yy+1),i)
                    B1 = self._read_data('b%i'%(subplot.xx+1),i)
                subplot.fieldlines = self.fieldlines(B1[subplot.arr_slice],B2[subplot.arr_slice],self.fieldseeds,dS=self.dS)
            self._draw_axes(subplot)
            self._draw_img(subplot)
            self.explots.append(subplot)
            subplot.axes.axis(subplot.extent)
        self.subplots = []
    
    def get_subplot_extent(self,subplot):
        #Get Required data
            arr_slice = subplot.arr_slice
            if not arr_slice == np.s_[...]:
                dims = np.array(list((type(s)!=int for s in arr_slice))).nonzero()[0]
                dims -= 2
                dims = np.abs(dims)
            else:
                dims = [0,1]
            subplot.yy = dims[1]
            subplot.xx = dims[0]
            subplot.xc = self.f.x[subplot.xx][arr_slice]
            subplot.yc = self.f.x[subplot.yy][arr_slice]
            #Define the extent of the grid
            extent = np.array((np.min(subplot.xc),
                                np.max(subplot.xc),
                                np.min(subplot.yc),
                                np.max(subplot.yc)))
            subplot.extent = extent
    
    def _draw_axes(self, subplot):
        subplot.axes = self.fig.add_subplot(self.nc,self.nr,subplot.subkey)
        subplot.axes.set_title(subplot.title)
        subplot.axes.axis(subplot.extent)
        
        mega = mpl.ticker.FuncFormatter(lambda x,pos:"$%1.1f$"%(x/1e6)) #TODO: this should be kwarg
        subplot.axes.xaxis.set_major_formatter(mega)
        subplot.axes.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
        subplot.axes.yaxis.set_major_formatter(mega)
        
    def _draw_img(self,subplot):
        """This handels the mpl plotting functions."""

        im = subplot.axes.imshow(subplot.data.T,origin='lower',
                                 extent=subplot.extent,norm=subplot.norm,
                                 cmap=subplot.cmap)

        scalar = mpl.ticker.ScalarFormatter(useMathText=False,useOffset=False)
        scalar.set_powerlimits((-3,3))
        subplot.axes.xaxis.set_major_formatter(scalar)
        subplot.axes.yaxis.set_major_formatter(scalar)
        subplot.x_offset = subplot.axes.xaxis.get_offset_text()
        subplot.x_offset.set_visible(False)
        subplot.y_offset  = subplot.axes.yaxis.get_offset_text()
        subplot.y_offset.set_visible(False)
        
        if subplot.labels and len(subplot.labels) == 2:
            subplot.axes.set_ylabel(subplot.labels[0])
            subplot.axes.set_xlabel(subplot.labels[1])
        elif subplot.labels:
            #TODO: Sort this out
            raise NotImplementedError("Label format not explicit... Sorry I am Being Lazy")
        
        if subplot.colorbar:
            divider = make_axes_locatable(subplot.axes)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(im,cax=cax)
            cbar.formatter = scalar
            if isinstance(subplot.norm,mpl.colors.LogNorm):
                cbar.locator = mpl.ticker.LogLocator()
            else:
                cbar.locator = mpl.ticker.MaxNLocator(9,prune='both')
            cbar.update_ticks()
            subplot.colorbar = cbar
        
        subplot.artists.append(im)
        
        if type(subplot.fieldlines) is np.ndarray:
            for line in subplot.fieldlines:
                line[0] = (line[0]/(len(subplot.xc)-2) * (subplot.xc.max() - subplot.xc.min()))
                line[1] = (line[1]/(len(subplot.yc)-2) * (subplot.yc.max() ))
                subplot.artists.append(subplot.axes.plot(line[0],line[1],color=self.line_colour,linewidth=1)[0])
        

    def seed_generator(self,ax,density,dS):
        """ Generate the fieldline seeds
        
            *ax* -Int (1,2,3,4) (top,bottom,left,right) - Side of grid to seed field lines
            *density* - Int (0-ng)- Field line density in grid points
            dS - Step size for fieldline calculation
        """
        #generte fixed seeds for all frames
        #Find number of points along chosen axis
        if ax == 1 or ax == 2:
            ng = self.f.header['nx'][0] -1 
        elif ax == 3 or ax == 4:
            ng = self.f.header['nx'][1] -1 
        else:
            raise Exception("Please specifiy axis as either 1,2,3,4")

        seeds = []
    
        if ax == 1:
            for each in np.arange(0,ng,density):
                seeds.append([each, float(ng)])
        elif ax == 2:
            for each in np.arange(0,ng,density):
                seeds.append([each, 0.])
        elif ax == 3:
            for each in np.arange(0,ng,density):
                seeds.append([0., each])
        elif ax == 4:
            for each in np.arange(0,ng,density):
                seeds.append([float(ng), each])
        else:
            raise Exception("Invalid axis specifier, also how did you get here")
        
        self.fieldseeds = seeds
        self.dS = dS
        return seeds
            
    def fieldlines(self,v1,v2,seeds,dS=1):
        """ Simple Euler fieldline integrator
        v1,v2 - y and x vectors
        seeds - array of coordinate (array indexes) pairs
        dS [optional] - step size
        """
        from scipy import ndimage
        
        field = []
        #Get extent of domain
        max1 = v1.shape[0]
        max2 = v2.shape[1]
        min2 = 0
        min1 = 0
        #integrate for each fieldline
        for seed in seeds:
            c1,c2 = seed
            out1 = [c1] #first point is seed
            out2 = [c2]
            cnt = 0
            
            while (c1 <= max1 and c1 >= min1) and (c2 <= max2 and c2 >= min2):
                #Interpolate the vector field to the coord
                coords = np.array([[c1],[c2]])
                v1i = ndimage.map_coordinates(v1,coords)[0]
                v2i = ndimage.map_coordinates(v2,coords)[0]
                vi = np.sqrt(v1i**2 + v2i**2)
                
                d1 = ( v1i * dS )/ vi
                d2 = ( v2i * dS )/ vi
                c1 -= d1 
                c2 -= d2
                
                out1.append(c1)
                out2.append(c2)
                
                cnt += 1
                if cnt > 500: # Maximum iteration limit
                    print "limit"
                    break
            out = np.zeros([len(out1),len(out2)])
            out[0] = out1
            out[1] = out2
            field.append(out)
        
        return np.array(field)

    def show(self):
        plt.show()

class SACanim(animation.TimedAnimation,SACplot):
    """This is a wrapper class for an animation routine."""
    def __init__(self, sacfile, figsize=(12,4), therange = None, sup_title='', interval = 200, blit = True):
        SACplot.__init__(self,sacfile,figsize=figsize,sup_title='')
        self.therange = therange or range(1, self.f.num_records,1) 
        
        #This needs to read in all data in "therange" and store it
        
        animation.TimedAnimation.__init__(self, self.fig, interval=interval, blit=blit)
        
    def animate(self):
        """Reads data, and plots requested subplots. Last routine to call."""
        self.fieldseq = []
        self.t = []

        keys = []
        for subplot in self.subplots:
            self.get_subplot_extent(subplot)
            self._draw_axes(subplot)
            subplot.anim_data = []
            keys.append(subplot.datakey)
            
        for i in self.therange:
            print "reading frame %i of %i"%(i,self.therange[-1])
            datas = self._read_data(keys,i)
            self.t.append(self.f.header['t'])
            for k,subplot in enumerate(self.subplots):
                #read all data sequence into subplot
                subplot.anim_data.append(datas[k][subplot.arr_slice])
                #Generate fieldlines for each frame
                #draw_img first frame
                if i == self.therange[0]:
                    subplot.data = subplot.anim_data[0]
                    self._draw_img(subplot)

            
            
    def _draw_frame(self, framedata):
        i = framedata
        for subplot in self.subplots:
            im = subplot.artists[0]
            im.set_data(subplot.anim_data[i].T)
            im.set_clim(np.min(subplot.anim_data[i]),np.max(subplot.anim_data[i]))
                
        
        
#        
        self.tbar[0].set_width(self.t[i])  
#        for k,art in enumerate(self.artists):
#            im = art[0]
#            im.set_data(self.data_w[i,k,:,:])
#            im.set_clim(np.min(self.data_w[i,k,:,:]),
#                              np.max(self.data_w[i,k,:,:]))
#
#            flines = art[1:]
#            if len(flines) > 0: #Catch fieldlines == False
#                for j,line in enumerate(self.fieldlines[i]):
#                    line = np.array(line)/(len(self.yc)/np.max(self.yc))
#                    flines[j].set_data(line)
                    
            
    def new_frame_seq(self):
        return iter(xrange(0,len(self.therange)))
        
    def _init_draw(self):
        pass

#if __name__ == "__main__":
#    import SACread
#    hslice = np.s_[:,:,6]
#    vslice = np.s_[:,62,:-4]
#    vlabels = ("$Z [Mm]$","$Y [Mm]$")
#    hlabels = ("$X [Mm]$", "$Y [Mm]$")
#    hfilename = "/home/stuart/iceberg_fastdata/3D_tube128_horiz_p30_C_A80000_driver30.h5"
#    hf = SACread.SAChdf5(hfilename)
#    Splt = SACanim(hf,sup_title="Hello",therange=range(1,20))
#    Splt.grid_spec(2,2)
#    Splt.add_img('rho',colorbar=True,arr_slice=vslice,norm='fixedcentre',cmap='jetblack')
#    Splt.add_img('v1',colorbar=True,title='$V_z [m/s]$',
#                 fieldlines=False,arr_slice=hslice,cmap='jetblack',norm='fixedcentre',labels=hlabels)
#    Splt.add_img('v2',colorbar=True,title='$V_x [m/s]$',
#                 fieldlines=False,arr_slice=vslice,cmap='jetblack',norm='fixedcentre',labels=vlabels)
#    Splt.add_img('v3',colorbar=True,title='$V_y [m/s]$',
#                 fieldlines=False,arr_slice=vslice,cmap='jetblack',norm='fixedcentre',labels=vlabels)    
#
##    Splt.plot(40)
##    Splt.add_tbar()
#    Splt.animate()
#    Splt.save("horiz_p30_C_A80000_anim.mp4", fps=5, codec='mpeg4', clear_temp=True,frame_prefix="/home/stuart/Documents/_tmp")
#    Splt.show()