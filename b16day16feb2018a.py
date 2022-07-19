from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/b16/data/2018/mt16638-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')

d(271879)

flipruns=range(271885,271945+7,7)
s1xvo=[];fracdiff=[]; err=[]
for run in flipruns:
    print d(run)
    s1xvo+=[mean(d.s1_s1vo)]
    fracdiff+=[mean(d.fracdiff)]
    err+=[std(d.fracdiff)]

figure();errorbar(s1xvo, fracdiff,err); grid(1);xlabel('s1vo'); ylabel('fracdiff')
#savefig('/scratch/plot1.pdf')





