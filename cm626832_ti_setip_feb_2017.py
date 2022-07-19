from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit
path='/dls/i16/data/2017/cm16772-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')


figure()
fluoscans=range(626832,626877)
for scan in fluoscans:
    d(scan)
    d.plot('energy2','sum',hold=1)

