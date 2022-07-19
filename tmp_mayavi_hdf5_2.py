###### get x,y,z data from hkl axes

import sys
sys.path
sys.path.append('/dls_sw/apps/scisoftpy/2.7')
sys.path.append('/dls_sw/i16/software/python')
#%matplotlib nbagg
from matplotlib.pyplot import *
matplotlib.verbose.set_level("helpful")
from numpy import *
import h5py
from mayavi import mlab
#mlab.init_notebook()

#f=h5py.File('/home/spc93/data/678450.h5','r') #read data
#f=h5py.File('/dls/science/users/spc93/processing/709614.h5','r') #read data
fdata=h5py.File('/dls/science/users/spc93/processing/709570.h5','r') #read data
fbg=h5py.File('/dls/science/users/spc93/processing/709582.h5','r') #read data
print fdata.keys()
#images_e1 = array(f.get('images_e1'))
#f.close()

ha=fdata.get('entry1/reciprocal_space/h-axis')
ka=fdata.get('entry1/reciprocal_space/k-axis')
la=fdata.get('entry1/reciprocal_space/l-axis')

v=fdata.get('entry1/reciprocal_space/volume')
vbg=fbg.get('entry1/reciprocal_space/volume')
vdat=array(v)
bgdat=array(vbg)

datnorm=vdat/vdat.max()
bgnorm=bgdat/bgdat.max()

#dat_minus_bg = datnorm - bgnorm
dat_minus_bg = datnorm

c=[2e-6]

yy, xx, zz = meshgrid(ka, ha, la)
#mlab.points3d(xx, yy, zz, dat_minus_bg)
mlab.contour3d(xx, yy, zz, dat_minus_bg, contours=c, opacity=.7)
mlab.axes(xlabel='h', ylabel='k', zlabel='l')
mlab.outline()
mlab.show()





