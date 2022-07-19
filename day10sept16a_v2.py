from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/i16/data/2016/cm14471-4/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')


figure()
d(612661) #energy scan - no sample - id calibrated at 10 keV
subplot(2,2,1); d.plot('energy2', 'ic1monitor', hold=1)
d.plot('energy2', 'diode', hold=1)

#compare ucal at 10 keV, 8.55 keV, 6.2 keV
subplot(2,2,3); d.plot('energy2', 'diode', hold=1)
d(612663); d.plot('energy2', 'diode', hold=1)
d(612665); d.plot('energy2', 'diode', hold=1)

figure()
d(612661)
d.plot('energy2', 'diode', hold=1)
d(612663); d.plot('energy2', 'diode', hold=1)
d(612665); d.plot('energy2', 'diode', hold=1)
savefig('/home/spc93/tmp/tmp.pdf')



figure()
#diamond powder scans
d(612675); diode_sample=d.diode; ic1_sample=d.ic1monitor
d(612676); diode_nosample=d.diode; ic1_nosample=d.ic1monitor
trans_rat=diode_sample/ic1_sample/diode_nosample*ic1_nosample
subplot(2,2,1); 
d.plot('energy2', trans_rat, hold=1)
subplot(2,2,3);
d(612675); d.plot('energy2', 'diode', hold=1)
d(612678); d.plot('energy2', 'diode', hold=1)


#figure()
#d(611898)#fluo scan
#subplot(2,2,1); d.plot('energy2','sum', hold=1)
#d(611934)#400 eta
#subplot(2,2,2); d.plot('eta','sum', hold=1)
#d(611940)#440 eta
#subplot(2,2,3); d.plot('eta','sum', hold=1)
#d(611962)#400 pol scan
#subplot(2,2,4); d.d(611898)plot('rotp','APD', hold=1)
#savefig('/home/spc93/tmp/tmp.pdf')

#figure()
#d(611964)#
#subplot(2,2,1); d.plot('pil100kthresh_keV',d.roi1_sum/max(d.roi1_sum), hold=1)
#d(611965)
#subplot(2,2,1); d.plot('pil100kthresh_keV',d.sum/max(d.sum), hold=1)
#d(611971)#vertical scattering 300
#subplot(2,2,2); d.plot('energy2',d.roi1_sum, hold=1)
#
#d(611994)#horizontal scattering 300 sig-pi psi~83
#subplot(2,2,3); d.plot('energy2',d.APD, hold=1)
#
#d(611997)#horizontal scattering 300 pilatus multiple scattering scan
#subplot(2,2,4); d.plot('psic','roi1_sum', hold=1)

#figure()
#for psiscan in range(612121,612197+1,2):
#    d(psiscan)
#    d.plot('psic','APD', hold=1)
#ylim(140000,250000)


#savefig('/home/i16user/tmp/tmp.pdf')