from dlstools import dataloader
import matplotlib.pyplot as plt
#import h5py
close('all')

datpath='/dls/i16/data/2017/cm16772-2/%i.dat'

d=dataloader.dlsloader(datpath)

d(641226)
d.plot('energy2','sum')

figure()
for scan in range(641250,641271):
    d(scan)
    plt.semilogy(d.l,d.APD, hold=1); 
plt.title('111 scans 6.53-6.55 keV'); grid(1); axis('tight')
#savefig('/home/spc93/tmp/plot.pdf')

figure()
d(641274)
cts=zeros(100)
for scan in range(641274,641293):
    d(scan)
    cts+=d.APD[0:100]
    plt.semilogy(d.l,d.APD, hold=1); 
plt.title('111 scans 6.5 keV'); grid(1); axis('tight')
figure(); plt.semilogy(d.l,cts, hold=1);  grid(1); axis('tight')
#savefig('/home/spc93/tmp/plot.pdf')

figure()
for scan in range(641326,641346+1):
    d(scan)
    plt.plot(d.l,d.APD, hold=1,label='%.4f keV' % d.en);
    #pause(1)
plt.title('111 scans over 11-1  6.53-6.55 keV'); grid(1); axis('tight')
legend(); plt.ylim([0,2000])
savefig('/home/spc93/tmp/plot.pdf')

figure()
for scan in range(641352,641363+1):
    d(scan)
    plt.plot(d.l,d.APD, hold=1,label='T= %.1f K' % mean(d.Ta));
plt.title('111 scans over 11-1  temp dependence'); grid(1); axis('tight')
legend(); plt.ylim([0,2000])
#savefig('/home/spc93/tmp/plot.pdf')

figure()
for scan in range(641401,641462+1,2):
    d(scan)
    plt.plot(d.l,d.APD, hold=1,label='%.4f keV' % d.en);
    pause(.5)
plt.title('hkl-tau 2 2 -2 along 111 vs energy'); grid(1); axis('tight')
legend(); plt.ylim([0,1000])
#savefig('/home/spc93/tmp/plot.pdf')

figure()
for scan in range(641401+1,641462+1,2):
    d(scan)
    plt.plot(d.l,d.APD, hold=1,label='%.4f keV' % d.en);
plt.title('hkl+tau 2 2 -2 along 111 vs energy'); grid(1); axis('tight')
legend(); plt.ylim([0,1000])
#savefig('/home/spc93/tmp/plot2.pdf')

figure()
for scan in range(641463,641486):
    d(scan)
    plt.plot(d.l,d.APD, hold=1,label='T= %.1f K' % mean(d.Ta));
plt.title('hkl-tau 2 2 -2 along 111 vs temperature going up'); grid(1); axis('tight')
legend(); plt.ylim([0,400])
savefig('/home/spc93/tmp/plot3.pdf')

figure()
for scan in range(641487,641500):
    d(scan)
    plt.plot(d.l,d.APD, hold=1,label='T= %.1f K' % mean(d.Ta));
plt.title('hkl-tau 2 2 -2 along 111 vs temperature going down'); grid(1); axis('tight')
legend(); plt.ylim([0,400])
savefig('/home/spc93/tmp/plot4.pdf')

