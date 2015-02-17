import h5py

from yt.frontends.gdf.data_structures import GDFDataset

from .fields import SACGDFFieldInfo

__all__ = ['SACGDFDataset']

class SACGDFDataset(GDFDataset):
    _field_info_class = SACGDFFieldInfo

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # If this is GDF, test it
        if super(SACGDFDataset, self)._is_valid(*args, **kwargs):
            try:
                fileh = h5py.File(args[0], 'r')
                data_software = fileh["gridded_data_format"].attrs['data_software']
                fileh.close()
            except:
                data_software = None
            return data_software in ('Sheffield Advanced Code', 'SAC')
        else:
            return False
