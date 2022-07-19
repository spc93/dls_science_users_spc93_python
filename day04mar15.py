from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit #need to fix - gca no longer in top namespace

close('all')

path='/dls/i16/data/2015/cm12169-1/'
d=dataloader.dlsloader(path+'%i.dat')
p=dataloader.tiffloader(d, lambda obj: path+obj.pilatus100k_path_template)


#figure(10);
#d(503926); subplot(2,2,1); d.plot('eta','sum','g', hold=1)
#d(503927); subplot(2,2,1); d.plot('eta','sum','r', hold=1)
#d(503928); subplot(2,2,1); d.plot('eta','sum','k', hold=1)
#d(503977); subplot(2,2,1); d.plot('Energy','sum', hold=1)

#figure(10);
#d(503979);
#for ii in range(1,4):
#    p(ii)
#    subplot(2,2,ii); pcolor(p.image_01); axis('tight'); title(p.file)
    #figure(); pcolor(p.image_01); axis('tight'); title(p.file)


#d(503984); d.plot('eta','sum')
#p(8); figure(); imshow(p.image_01); axis('tight'); title(p.file)


#d(503995)
#figure(10); 
#subplot(2,2,1); d.plot('Energy','sum', hold=1)
#d(504009)
#subplot(2,2,1); d.plot('Energy','sum', hold=1)

#d(504010); p(240)
#p.plot()

#figure(10); 
#subplot(2,2,1); d.plot('Energy','sum', hold=1)

#ss=d.findscans(range(504037,505117),'self.Energy>7.2999')
#d.sequence(range(504038,504237+1)) #eta vs energy
#figure(); gca(projection='3d').plot_wireframe(d.s.Energy, d.s.eta, d.s.sum, label='GGG Gd L3 008')

area=[]; centre=[]; area_err=[]; centre_err=[]; width=[]; width_err=[]; en=[];
for n in range(504038,504237):
    d(n)
    figure(10)
    d.plot('eta','sum',hold=0)
    fit(pv_c)
    pause(0.1)
    area+=[pv_c.p[0]]; area_err+=[pv_c.err[0]]
    centre+=[pv_c.p[1]-mean(d.eta)]; centre_err+=[pv_c.err[1]]
    width+=[pv_c.p[2]]; width_err+=[pv_c.err[2]];
    en+=[d.Energy]
area=array(area); centre=array(centre); area_err=array(area_err); centre_err=array(centre_err); en=array(en)
subplot(2,2,1); errorbar(en, centre, yerr=centre_err); xlabel('Energy (keV)'); ylabel('Theta (deg)'); axis('tight');
subplot(2,2,2); errorbar(en, area, yerr=area_err); xlabel('Energy (keV)'); ylabel('Area'); axis('tight');
subplot(2,2,3); errorbar(en, width, yerr=width_err); xlabel('Energy (keV)'); ylabel('Width (deg)'); axis('tight');


# stability: range(504238,505118)
area=[]; centre=[]; area_err=[]; centre_err=[]; width=[]; width_err=[]; hours=[]; scn=[];
for n in range(504239,505117):
    d(n)
    figure(11)
    d.plot('eta','sum',hold=0)
    fit(pv_c)
    pause(0.1)
    area+=[pv_c.p[0]]; area_err+=[pv_c.err[0]]
    centre+=[pv_c.p[1]-mean(d.eta)]; centre_err+=[pv_c.err[1]]
    width+=[pv_c.p[2]]; width_err+=[pv_c.err[2]];
    scn+=[n]
area=array(area); centre=array(centre); area_err=array(area_err); centre_err=array(centre_err); en=array(en); scn=array(scn);
subplot(2,2,1); errorbar(scn, centre, yerr=centre_err); xlabel('scan number'); ylabel('Theta (deg)'); axis('tight');
subplot(2,2,2); errorbar(scn, area, yerr=area_err); xlabel('scan number'); ylabel('Area'); axis('tight');
subplot(2,2,3); errorbar(scn, width, yerr=width_err); xlabel('scan number'); ylabel('Width (deg)'); axis('tight');


#savefig('/home/i16user/tmp/tmp.pdf')



