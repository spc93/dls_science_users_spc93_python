from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit #need to fix - gca no longer in top namespace

close('all')

path='/dls/i16/data/2016/cm14471-4/'
d=dataloader.dlsloader(path+'%i.dat')




d(613672); subplot(2,2,1); d.plot('delta','APD', hold=1)
d(613673); subplot(2,2,1); d.plot('delta','APD', hold=1)