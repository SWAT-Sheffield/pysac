import mhs_atmosphere as mhs
import yt
ds=yt.load(mhs.filename)
cg=ds.index.grids[0]
print 'rho.max cgs =',cg['density_bg'].max()
print 'rho.max SI =',mhs.rho.max()
