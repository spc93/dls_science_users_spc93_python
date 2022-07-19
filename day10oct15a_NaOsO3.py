from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit
path='/dls/i16/data/2015/cm12169-4/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')



d(536211)#400 eta
subplot(2,2,1); d.plot('eta','sum', hold=1)
d(536219)#440 eta
subplot(2,2,2); d.plot('eta','sum', hold=1)
d(536223)#L3 fluo
subplot(2,2,3); d.plot('energy2','sum', hold=1)
figure()
d(536226)#300 ATS psi 90
subplot(2,2,1); d.plot('eta','roi1_sum', hold=1)
d(536227)#500 ATS psi 90
subplot(2,2,2); d.plot('eta','roi1_sum', hold=1)
d(536232)#pil100kthresh (300)

#pil100kgain(2) (high) - threshold max ~ 10
subplot(2,2,3); d.plot('pil100kthresh_keV',d.sum/max(d.sum), hold=1)
d(536235)#pil100kthresh (400, atten 255)
subplot(2,2,3); d.plot('pil100kthresh_keV',d.sum/max(d.sum), hold=1)
d(536236)#pil100kthresh (400, atten 255, off eta to reduce count rate)
subplot(2,2,3); d.plot('pil100kthresh_keV',d.sum/max(d.sum), hold=1)
d(536238)#pil100kthresh (400, atten 255, off eta more to reduce count rate and use scroi)
subplot(2,2,3); d.plot('pil100kthresh_keV',d.roi1_sum/max(d.roi1_sum), hold=1); axis('tight')

#pil100kgain(1) (med)
d(536239)#pil100kthresh (400, atten 255)
subplot(2,2,4); d.plot('pil100kthresh_keV',d.sum/max(d.sum), hold=1)
d(536240)#pil100kthresh (fluo, atten 0)
subplot(2,2,4); d.plot('pil100kthresh_keV',d.sum/max(d.sum), hold=1)
axis('tight')

figure()

d(536241)#300 ATS psi 90 with thresh=10 and background roi
subplot(2,2,1); d.plot('eta','roi1_sum', hold=1)
subplot(2,2,1); d.plot('eta','roi1_bg_sum', hold=1)

figure()
for i in range(536252,536259):
    d(i); d.plot('psic','roi1_sum', hold=1)

figure()
for i in range(536264,536277,3): #300 short energy scans
    d(i); d.plot('energy2','roi1_sum', hold=1)
    d.plot('energy2','roi1_bg_sum','r', hold=1)

figure(99); subplot(2,1,1)
for i in range(536282,536297+1,3): #300 long energy scans
    d(i); d.plot('energy2','roi1_sum', hold=1)
    d.plot('energy2','roi1_bg_sum','r', hold=1)
axis('tight'); ylim([0,120000]); xlabel('Energy (keV)')

figure()
psi=[]; pkarea=[]
for scan in range(536298,536448+1): #300 azimuthal scans
    d(scan)
    psi+=[d.psi+360*(d.psi<0)]#add 360 deg if -ve
    cts=d.roi1_sum-d.roi1_bg_sum
    d.plot('eta',cts, hold=1 ); ylim([0,150000])
    [centre, width, sum, height, area, m, c, bg]=peak(d.eta,cts ,nbgpts=5)
    #pause(.1)
    pkarea+=[area]
psi=array(psi); pkarea=array(pkarea);
figure(99); subplot(2,1,2); 
plot(psi, pkarea,'go'); ylim([-100,10000]); grid(1); xlabel('$\psi$ (ref = 001)'); ylabel('Peak area')

savefig('/home/spc93/tmp/tmp.pdf')