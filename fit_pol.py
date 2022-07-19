from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit
from polarization_analyser_formula import PA

close('all')


energy_keV=10   #use if no th_pol_deg given
d_spacing=4     #use if no th_pol_deg given
th_pol_deg=45   #angle if known, else None
if th_pol_deg==None:
    th=arcsin(6.2/d_spacing/energy_keV)
    th_pol_deg=180/pi*th

d=dataloader.dlsloader('/dls/i16/data/2015/cm12169-4/%i.dat')

PA.params['theta_pol'].vary=True
PA.params['StokesP3'].vary=False
PA.params['StokesP1'].vary=True
d(535440)
plot(d.stoke,d.APD)
PA.fit(d.stoke, d.APD, pin=[35,1e6,0,0.99])
print PA
hold(1)
plot(d.stoke,PA(d.stoke))



PA.params['theta_pol'].vary=False
PA.params['StokesP3'].vary=True
PA.params['StokesP1'].vary=True

pppscan1=range(535463, 535503+1)
pppscan2=range(535504, 535885)


figure(); hold(1)
pp_offet=[]; P1=[]; P2=[]; P3=[];
for pp in pppscan1:
    d(pp)
    pp_offet+=[mean(d.ppp_offset)]
    plot(d.stoke,d.APD)
    PA.fit(d.stoke, d.APD, pin=[35.5,1,1,0])
    print PA
    hold(1); plot(d.stoke, PA(d.stoke)); axis('tight')
    P1+=[PA.p[2]]; P2+=[sqrt(1-PA.p[2]**2-PA.p[3]**2)]; P3+=[PA.p[3]]; 
    pause(.5); print mean(d.ppp_offset) 
figure(); plot(pp_offet, P1, 'ro', hold=1, label='P1',); plot(pp_offet, P2, 'go', hold=1, label='Not linear'); plot(pp_offet, P3, 'ko', hold=1, label='P3',); axis('tight'); grid(1); legend()
xlabel('ppp offset (degrees)'); ylabel('Stokes parameters')

savefig('/home/spc93/tmp/PA_pppa111_11p564keV_short.pdf')
#savefig('/home/spc93/tmp/PA_pppa111_11p564keV.pdf')
