from dlstools.dataloader import pdnx
from dlstools import dataloader
import pandas as pd

close('all')

p='/dls/i16/data/2017/mt16226-1/%i.nxs'


dataloader._scandata_field='/entry1/pil100k'
n=pdnx(p % 637131)
print n.scan
n.plot()
plot(n.eta, n.nx.entry1.roi1.roi1_sum.nxdata)


from dlstools import dataloader
path='/dls/i16/data/2017/mt16226-1/'
d=dataloader.dlsloader(path+'%i.dat')
p=dataloader.piloader(d, lambda obj: path+obj.andor_path_template)


#for scan in range(637200,637299+1):
#    try:
#        print d(scan); print d.Tset
#    except:
#        print "=== error reading file"



#d(637247); p(1);
#roi_ij=[150,300,300,600] 
#sum_im=zeros(p.image_01[roi_ij[0]:roi_ij[1], roi_ij[2]:roi_ij[3]].shape)
#for i in d.path:
#    p(i)
#    sum_im+=p.image_01[roi_ij[0]:roi_ij[1], roi_ij[2]:roi_ij[3]]
#figure(); imshow(sum_im); axis('tight'); title(p.file)

scansup=range(637247+1,637272+2,2)
scansdown=range(637273+1,637298+2,2)
d(637247); p(1);
roi_ij=[150,300,300,600] 

figure()
for scan in scansup+scansdown:
    d(scan)
    sum_im=zeros(p.image_01[roi_ij[0]:roi_ij[1], roi_ij[2]:roi_ij[3]].shape)
    for i in d.path:
        p(i); sum_im+=p.image_01[roi_ij[0]:roi_ij[1], roi_ij[2]:roi_ij[3]]
    #figure(); 
    imshow(sum_im, clim=(1.5e5,4e5)); axis('tight'); title('T = '+str(d.Tset)+' K'); 
    #colorbar()
    savefig('/home/spc93/tmp/acme%i.jpg' % scan)
#    figure(); imshow(sum_im); axis('tight'); title('T = '+str(d.Tset)+' K')

