from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit
path='/dls/i16/data/2014/mt10132-1/'

close('all')

#create .dat file loader object
d=dataloader.dlsloader(path+'%i.dat')
#create .tiff file loader object which refers to current .dat file
p=dataloader.tiffloader(d, lambda obj: path+obj.pilatus100k_path_template)

#close('all') to close all plot windows
#load a .dat file...
print d(477007) #using print is optional and displays scan info from within script
#load the second image from it...
print p(20)
figure(); imshow(p.image_01); axis('tight'); title(p.file)

#load another .dat file...
d(477010)
#and now load the next one...
d.inc()
#quick way to plot with axis labels (see help: d.plot?)
d.plot('eta','sum')

#fit pseudo Voigt
fit(pv_c)

#example plot simple 3d scan
#figure().add_subplot(111, projection='3d').plot_wireframe(d.rotp, d.thp, d.APD)

#examples of 3d plots from sequence of scans
#d.sequence(range(541522,541542),['thp','stoke','APD'])
#figure().add_subplot(111, projection='3d').plot_wireframe(d.s.stoke, d.s.thp, d.s.APD); axis('tight')
#figure().add_subplot(111, projection='3d').plot_surface(d.s.stoke, d.s.thp, d.s.APD, cmap='autumn', cstride=1, rstride=1); axis('tight')
