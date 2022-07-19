from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit
close('all')

path='/dls/i16/data/2010/mt1690-1/'
d=dataloader.dlsloader(path+'%i.dat')

print d(127750)
figure(); plot(d.Energy,d.sum/max(d.sum)); xlabel('Energy (keV)'); ylabel('Fluo(norm)'); title('#127750 fluo 100 (010) wedge'); axis('tight'); grid(1)
xlim([4.96,4.99]); grid(1);

print d(127577); figure(); plot(d.Energy,-log(d.APD))
hold(1); d(127580); plot(d.Energy,-log(d.APD),'r')
hold(1); d(127583); plot(d.Energy,-log(d.APD),'g')
hold(1); d(127586); plot(d.Energy,-log(d.APD),'c')
title('Borrmann 200 (010) wedge 2.5 mm from top')
axis('tight'); xlim([4.96,4.99]); grid(1);
print 'pre-edge peak ratios look ~ same for Borrmann - all enhanced by same factor'
print 'not sure what this means. 200 borrmann looks very odd anyway' 
#p=dataloader.tiffloader(d, lambda obj: path+obj.pilatus100k_path_template)

path='/dls/dmf/i16/data/2007/mt0'
d=dataloader.dlsloader(path+'%i.dat')
d(53629)



#print d(477007)
#print p(2)
#figure(); pcolor(p.image0); axis('tight'); title(p.file)
#print d(477000)
#figure(); plot(d.eta, d.sum); title(d.file+'\n'+d.cmd); axis('tight'); grid(1)
#fit(pv_c)



#d(482940)
#subplot(2,2,1); plot(d.Energy, d.sum); grid(1); xlabel('Energy'); ylabel('sum'); title(d.file)
#savefig('/home/i16user/tmp/tmp.pdf')