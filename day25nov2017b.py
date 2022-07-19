#TiO2 dichroism to check angles

from dlstools import dataloader
close('all')
path='/dls/i16/data/2017/cm16772-4/'
d=dataloader.dlsloader(path+'%i.dat')


scans=range(664361,664432+1)
#scans=range(664361,664365)

kphi=[]
for scan in scans:
    d(scan)
    kphi+=[d.kphi[0]]
    d.plot('energy2','sum',hold=True)

savefig('/home/spc93/tmp/tio2_spec.pdf')

d.sequence(scans)
figure().add_subplot(111, projection='3d').plot_wireframe(d.s.kphi, d.s.energy2, d.s.sum); axis('tight')
savefig('/home/spc93/tmp/tio2_wire.pdf')
