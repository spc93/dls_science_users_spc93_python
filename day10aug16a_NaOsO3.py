from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/i16/data/2016/cm14471-4/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')

#figure()
#d(611898)#fluo scan
#subplot(2,2,1); d.plot('energy2','sum', hold=1)
#d(611934)#400 eta
#subplot(2,2,2); d.plot('eta','sum', hold=1)
#d(611940)#440 eta
#subplot(2,2,3); d.plot('eta','sum', hold=1)
#d(611962)#400 pol scan
#subplot(2,2,4); d.plot('rotp','APD', hold=1)
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



#last=612454
mu1scans=array([612210]+range(612228,612311+1,6)+range(612317,612454+1,6))
#exception raised #612316
gamscans=mu1scans+1
chiscans=mu1scans+2
mu2scans=mu1scans+3
e90scans=mu1scans+4
e0scans=mu1scans+5

figure()
subplot(2,2,1);
for scan in mu1scans:
    d(scan);    d.plot('mu','APD', hold=1)
#subplot(2,2,2);
#for scan in gamscans:
#    d(scan);    d.plot('gam','APD', hold=1)
#subplot(2,2,3);
#for scan in chiscans:
#    d(scan);    d.plot('chi','APD', hold=1)
subplot(2,2,4);
for scan in mu2scans:
    d(scan);    d.plot('mu','APD', hold=1)

figure()
subplot(2,2,1);
for scan in e90scans:
    d(scan);    d.plot('energy2','APD', hold=1)
subplot(2,2,2);
for scan in e0scans:
    d(scan);    d.plot('energy2','APD', hold=1)

#subplot(2,2,1); d.plot('pil100kthresh_keV',d.roi1_sum/max(d.roi1_sum), hold=1)
#d(611965)

savefig('/home/spc93/tmp/tmp.pdf')


fig = plt.figure()
ax = fig.gca(projection='3d')
d.sequence(mu2scans,'psi')
ax.plot_wireframe(d.s.mu, d.s.psi, d.s.APD, cstride=0)
xlabel('theta'); ylabel('Azimuth');
title('300 theta scans pi-sigma')
savefig('/home/spc93/tmp/tmp_th.pdf')

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#d.sequence(gamscans,'psi')
#ax.plot_wireframe(d.s.gam, d.s.psi, d.s.APD, cstride=0)

fig = plt.figure()
ax = fig.gca(projection='3d')
d.sequence(e90scans,'psi')
ax.plot_wireframe(d.s.energy2, d.s.psi, d.s.APD, cstride=0)
axis('tight'); ax.set_zlim3d(0,20000)
xlabel('Energy(keV)'); ylabel('Azimuth');
title('300 Energy scans pi-sigma')
savefig('/home/spc93/tmp/tmp_e90.pdf')

fig = plt.figure()
ax = fig.gca(projection='3d')
d.sequence(e0scans,'psi')
ax.plot_wireframe(d.s.energy2, d.s.psi, d.s.APD, cstride=0)
axis('tight'); ax.set_zlim3d(0,2000)
xlabel('Energy(keV)'); ylabel('Azimuth');
title('300 Energy scans pi-pi')
savefig('/home/spc93/tmp/tmp_e0.pdf')