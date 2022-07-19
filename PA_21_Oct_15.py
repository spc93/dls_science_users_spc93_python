from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit

close('all')

d=dataloader.dlsloader('/dls/i16/data/2015/cm12169-4/%i.dat')

#d(541516)# single 2d scan
#figure().add_subplot(111, projection='3d').plot_wireframe(d.rotp, d.thp, d.APD)


d.sequence(range(541522,541542),['thp','stoke','APD'])
figure().add_subplot(111, projection='3d').plot_surface(d.s.stoke, d.s.thp, d.s.APD, cmap='autumn', cstride=1, rstride=1); axis('tight')

d(541543); d.plot('rotp','APD')