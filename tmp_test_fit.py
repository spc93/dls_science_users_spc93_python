datadir='/dls/b16/data/2018/mt17690-1/'

import sys
sys.path
sys.path.append('/dls_sw/apps/scisoftpy/2.7')
sys.path.append('/dls_sw/i16/software/python')
from matplotlib.pyplot import *
#%matplotlib nbagg
from numpy import *
from dlstools import dataloader
from dlstools.dataloader import vec2mat


from dlstools.quickfit import *
from dlstools.dirty_fit import fit

d=dataloader.dlsloader(datadir+'%i.dat')
p=dataloader.tiffloader(d, lambda obj: datadir+obj.pilatus2m_path_template)


d(265322)
d.plot('mtYaw', 'sum')
#fit(gauss_c)
