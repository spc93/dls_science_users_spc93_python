from dlstools import dataloader
import matplotlib.pyplot as plt
#import h5py
close('all')

datpath='/dls/i16/data/2017/cm16772-2/%i.dat'

d=dataloader.dlsloader(datpath)

figure()
for scan in range(641110, 641115+1):
    d(scan)
    d.plot('en','ic1', hold=1)
savefig('/home/spc93/tmp/offset.pdf')


figure()
scans=[641118, 641119, 641121, 641122, 641123, 641124, 641125]
for scan in scans:
    d(scan)
    d.plot('en',d.ic1/d.rc, hold=1)
savefig('/home/spc93/tmp/taper.pdf')
