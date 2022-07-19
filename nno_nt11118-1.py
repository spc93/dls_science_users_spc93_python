from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
close('all')
path='/dls/i16/data/2014/nt11118-1/'
#create .dat file loader object
d=dataloader.dlsloader(path+'%i.dat')
#create .tiff file loader object which refers to current .dat file
p=dataloader.tiffloader(d, lambda obj: path+obj.pilatus100k_path_template)


d(474791)
d.plot('eta','sum')
d.plot('path','sum')
p(19)
figure(); pcolor(p.image0); axis('tight'); title(p.file); clim([0,10])

p(29)
figure(); pcolor(p.image0); axis('tight'); title(p.file); clim([0,1000])