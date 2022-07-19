from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit #need to fix - gca no longer in top namespace




#All coincidence events of ToA time difference {-150*2ns & 150*2ns}.
#Collection Time : 62 min
d=np.loadtxt('/dls_sw/i16/software/python/userscripts/steve/All Events Time Difference Histogram.txt').transpose()

figure()
subplot(2,2,1)
plot(d[0]*2,d[1]); title('TimePix coincidence spectrum');
xlabel('Time (ns)'); ylabel('Counts')
savefig('/home/spc93/tmp/timepix1.pdf')