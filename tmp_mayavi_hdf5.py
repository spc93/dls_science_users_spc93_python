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
f=h5py.File('/dls/science/users/spc93/processing/709576.h5','r') #read data

print f.keys()
#images_e1 = array(f.get('images_e1'))
#f.close()

ha=f.get('entry1/reciprocal_space/h-axis')
ka=f.get('entry1/reciprocal_space/k-axis')
la=f.get('entry1/reciprocal_space/l-axis')

v=f.get('entry1/reciprocal_space/volume')

vdat=array(v)

#works fine
#from numpy import pi, sin, cos, mgrid
#dphi, dtheta = pi/250.0, pi/250.0
#[phi,theta] = mgrid[0:pi+dphi*1.5:dphi,0:2*pi+dtheta*1.5:dtheta]
#m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
#r = sin(m0*phi)**m1 + cos(m2*phi)**m3 + sin(m4*theta)**m5 + cos(m6*theta)**m7
#x = r*sin(phi)*cos(theta)
#y = r*cos(phi)
#z = r*sin(phi)*sin(theta)
# View it.
#from mayavi import mlab
#s = mlab.mesh(x, y, z)
#mlab.show()

#mlab.contour3d(log(vdat+1), contours=4, transparent=True)
#c=[vdat.max()/2, vdat.max()/8, vdat.max()/32, vdat.max()/128, ]
#c=[vdat.max()/2, vdat.max()/32 ]
c=[vdat.max()/2, vdat.max()/200000 ]

#mlab.contour3d(vdat, contours=c, opacity=0.4)
#mlab.axes()
#mlab.show()

yy, xx, zz = meshgrid(ka, ha, la)
mlab.contour3d(xx, yy, zz, vdat, contours=c, opacity=.7)
mlab.axes(xlabel='h', ylabel='k', zlabel='l')
mlab.outline()
mlab.show()





