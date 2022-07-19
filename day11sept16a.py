from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/i16/data/2016/cm14471-4/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')


figure()

subplot(2,2,1); 
d(612708); d.plot('l', 'roi2_sum', hold=1)

subplot(2,2,3); 
d(612711); d.plot('k', 'roi2_sum', hold=1)


figure()
subplot(2,2,1); 
d(612714); d.plot('h', 'roi2_sum', hold=1)

figure()
subplot(2,2,1); 
for scan in range(612719,612753-2,2):
    d(scan); d.plot('idgap', 'ic1monitor', 'g',hold=1)
    d.inc(); d.plot('idgap', 'ic1monitor', 'r',hold=1)

subplot(2,2,3); #repeat and zoom in
for scan in range(612719,612753-2,2):
    d(scan); d.plot('idgap', 'ic1monitor', 'g',hold=1)
    d.inc(); d.plot('idgap', 'ic1monitor', 'r',hold=1)


savefig('/home/i16user/tmp/tmp.pdf')