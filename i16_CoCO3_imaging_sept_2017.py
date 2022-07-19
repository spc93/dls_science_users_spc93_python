from dlstools import dataloader
from quickfit import *
path='/dls/i16/data/2017/mt16227-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')

d(652877)

Ta=[]; Tb=[]; Dsum=[]
#for n in range(652877,653025,2):
for n in range(653025,653172,2):
    d(n)
    [centre, width, sum, height, area, m, c, bg]=peak(d.eta,d.sum,10)
    Ta+=[d.Ta]; Tb+=[d.Tb]; Dsum+=[sum] 
plt.plot(Ta, Dsum,'og',hold=1)

Ta=[]; Tb=[]; Dsum=[]
#for n in range(652877+1,653025,2):
for n in range(653025+1,653172,2):
    d(n)
    [centre, width, sum, height, area, m, c, bg]=peak(d.eta,d.sum,10)
    Ta+=[d.Ta]; Tb+=[d.Tb]; Dsum+=[sum] 
plt.plot(Ta, Dsum,'or',hold=1)
grid(1)

savefig('/dls/science/users/spc93/tmp/tmp.pdf')


'''

samples=['diamond_6um','diamond_3um','silicon_nbs640']
for sample in samples:
    figure()
    for scan in range(613679, 613993+1):
        d(scan)
        if sample in d.note:
            if 'Bragg edge' in d.note:
                rat_edge=d.diode*1.0/d.IC2
            if 'no sample' in d.note:
                rat_null=d.diode*1.0/d.IC2
                double_ratio=rat_edge/rat_null
                subplot(2,2,1); plot(d.energy2, double_ratio, hold=1); axis('tight'); title('Sample = '+sample)
                subplot(2,2,2); plot(d.energy2-mean(d.energy2), double_ratio/mean(double_ratio), hold=1); axis('tight'); title('Sample = '+sample)
            if 'Bragg peak' in d.note:
                subplot(2,2,3); plot(d.energy2, d.APD, hold=1); axis('tight'); title('Sample = '+sample)
                subplot(2,2,4); plot(d.energy2-mean(d.energy2), d.APD*1.0/max(d.APD), hold=1); axis('tight'); title('Sample = '+sample)
    savefig('/home/spc93/tmp/%s.pdf' % sample)



#plot(d[0]*2,d[1]); title('TimePix coincidence spectrum');
#xlabel('Time (ns)'); ylabel('Counts')
#savefig('/home/spc93/tmp/tmp.pdf')

'''