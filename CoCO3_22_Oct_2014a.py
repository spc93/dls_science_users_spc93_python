#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as fitplt
from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit

path='/dls/i16/data/2014/cm4968-5/'
d=dataloader.dlsloader(path+'%i.dat')
p=dataloader.tiffloader(d, lambda obj: path+obj.pilatus100k_path_template)


#print d(477007)
#print p(2)
#figure(); pcolor(p.image0); axis('tight'); title(p.file)
#print d(477000)
#figure(); plot(d.eta, d.sum); title(d.file+'\n'+d.cmd); axis('tight'); grid(1)
#fit(pv_c)



d(482940)
subplot(2,2,1); plot(d.Energy, d.sum); grid(1); xlabel('Energy'); ylabel('sum'); title(d.file)
#savefig('/home/i16user/tmp/tmp.pdf')