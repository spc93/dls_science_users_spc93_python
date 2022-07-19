#focus and flux vs s1 slit size

from dlstools import dataloader
from matplotlib.pyplot import subplots_adjust
from pyqtgraph.examples.customGraphicsItem import plt
close('all')
path='/dls/i16/data/2017/cm16772-4/'
d=dataloader.dlsloader(path+'%i.dat')

M=25.5  #mono dist
M1=30.0   #mirror1
M2=31.8   #mirror2
F=50    #focus
#assume focus changes by demag of movement of virtual source
#Dt=delta theta (mono) (vert)
#Dt=delta theta (mono)/sin(theta) (horiz) (?)
#Dfy=M*(1-M1/F)*Dt #vert
#Dfx=M*(1-M2/F)*Dt/sin(theta) #horiz

#calculate Dt vs s1 xgap for s1ygap=1

print '== si thermal expansion crosses zero at T=120K'
print '=== check thermal conductivity for natural and single isotope'
print '=== operate at 120K??'

figure()
d(664449)
idx=find(d.s1ygap==1)
Dt=d.fwhmy_mm/M/(1-M1/F)*1000 #urad
subplot(2,2,1); plot(d.s1xgap[idx], d.sum[idx]); xlabel('xgap'); ylabel('sum'); grid(); title('harmonic=3')
subplot(2,2,2); plot(d.s1xgap[idx], d.T1dcmSi111[idx]); xlabel('xgap'); ylabel('T1dcm'); grid()
subplot(2,2,3); plot(d.s1xgap[idx], d.fwhmy_mm[idx]); xlabel('xgap'); ylabel('fwhm_y'); grid()
subplot(2,2,4); plot(d.s1xgap[idx], Dt[idx]); xlabel('xgap'); ylabel('slope_v_urad'); grid()
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1xgap_uharmonic3.pdf')

figure()
d(664452)
idx=find(d.s1ygap==1)
Dt=d.fwhmy_mm/M/(1-M1/F)*1000 #urad
subplot(2,2,1); plot(d.s1xgap[idx], d.sum[idx]); xlabel('xgap'); ylabel('sum'); grid(); title('harmonic=5')
subplot(2,2,2); plot(d.s1xgap[idx], d.T1dcmSi111[idx]); xlabel('xgap'); ylabel('T1dcm'); grid()
subplot(2,2,3); plot(d.s1xgap[idx], d.fwhmy_mm[idx]); xlabel('xgap'); ylabel('fwhm_y'); grid()
subplot(2,2,4); plot(d.s1xgap[idx], Dt[idx]); xlabel('xgap'); ylabel('slope_v_urad'); grid()
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1xgap_uharmonic5.pdf')

figure()
d(664455)
idx=find(d.s1ygap==1)
Dt=d.fwhmy_mm/M/(1-M1/F)*1000 #urad
subplot(2,2,1); plot(d.s1xgap[idx], d.sum[idx]); xlabel('xgap'); ylabel('sum'); grid(); title('harmonic=7')
subplot(2,2,2); plot(d.s1xgap[idx], d.T1dcmSi111[idx]); xlabel('xgap'); ylabel('T1dcm'); grid()
subplot(2,2,3); plot(d.s1xgap[idx], d.fwhmy_mm[idx]); xlabel('xgap'); ylabel('fwhm_y'); grid()
subplot(2,2,4); plot(d.s1xgap[idx], Dt[idx]); xlabel('xgap'); ylabel('slope_v_urad'); grid()
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1xgap_uharmonic7.pdf')

figure()
d(664458)
idx=find(d.s1ygap==1)
Dt=d.fwhmy_mm/M/(1-M1/F)*1000 #urad
subplot(2,2,1); plot(d.s1xgap[idx], d.sum[idx]); xlabel('xgap'); ylabel('sum'); grid(); title('harmonic=9')
subplot(2,2,2); plot(d.s1xgap[idx], d.T1dcmSi111[idx]); xlabel('xgap'); ylabel('T1dcm'); grid()
subplot(2,2,3); plot(d.s1xgap[idx], d.fwhmy_mm[idx]); xlabel('xgap'); ylabel('fwhm_y'); grid()
subplot(2,2,4); plot(d.s1xgap[idx], Dt[idx]); xlabel('xgap'); ylabel('slope_v_urad'); grid()
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1xgap_uharmonic9.pdf')








'''
d(664449);figure(); 
subplot(2,2,1); d.plot('T1dcmSi111','fwhmx_mm','.g',hold=1); ylim([.1,.2]); grid(1)
subplot(2,2,2); d.plot('T1dcmSi111','fwhmy_mm','.g',hold=1); ylim([0.02,.05]); grid(1)
d(664452);
subplot(2,2,3); d.plot('T1dcmSi111','fwhmx_mm','.r',hold=1); ylim([.1,.2]); grid(1)
subplot(2,2,4); d.plot('T1dcmSi111','fwhmy_mm','.r',hold=1); ylim([0.02,.05]); grid(1)
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1size_uharmonic35_t1dcm.pdf')

d(664455);figure(); 
subplot(2,2,1); d.plot('T1dcmSi111','fwhmx_mm','.g',hold=1); ylim([.1,.22]); grid(1)
subplot(2,2,2); d.plot('T1dcmSi111','fwhmy_mm','.g',hold=1); ylim([0.02,.1]); grid(1)
d(664458);
subplot(2,2,3); d.plot('T1dcmSi111','fwhmx_mm','.r',hold=1); ylim([.1,.22]); grid(1)
subplot(2,2,4); d.plot('T1dcmSi111','fwhmy_mm','.r',hold=1); ylim([0.02,.1]); grid(1)
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1size_uharmonic79_t1dcm.pdf')


d(664449);figure(); 
d.dict['beam_area']=d.s1xgap*d.s1ygap;
subplot(2,2,1); d.plot('sum','fwhmx_mm','.g',hold=1); ylim([.1,.2]); grid(1)
subplot(2,2,2); d.plot('sum','fwhmy_mm','.g',hold=1); ylim([0.02,.05]); grid(1)
subplot(2,2,3); d.plot('beam_area','fwhmx_mm','.r',hold=1); ylim([.1,.2]); grid(1)
subplot(2,2,4); d.plot('beam_area','fwhmy_mm','.r',hold=1); ylim([0.02,.05]); grid(1)
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1size_uharmonic3.pdf')


d(664452);figure(); 
d.dict['beam_area']=d.s1xgap*d.s1ygap;
subplot(2,2,1); d.plot('sum','fwhmx_mm','.g',hold=1); ylim([.1,.22]); grid(1)
subplot(2,2,2); d.plot('sum','fwhmy_mm','.g',hold=1); ylim([0.02,.08]); grid(1)
subplot(2,2,3); d.plot('beam_area','fwhmx_mm','.r',hold=1); ylim([.1,.22]); grid(1)
subplot(2,2,4); d.plot('beam_area','fwhmy_mm','.r',hold=1); ylim([0.02,.08]); grid(1)
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1size_uharmonic5.pdf')


d(664455);figure(); 
d.dict['beam_area']=d.s1xgap*d.s1ygap;
subplot(2,2,1); d.plot('sum','fwhmx_mm','.g',hold=1); ylim([.1,.22]); grid(1)
subplot(2,2,2); d.plot('sum','fwhmy_mm','.g',hold=1); ylim([0.02,.1]); grid(1)
subplot(2,2,3); d.plot('beam_area','fwhmx_mm','.r',hold=1); ylim([.1,.22]); grid(1)
subplot(2,2,4); d.plot('beam_area','fwhmy_mm','.r',hold=1); ylim([0.02,.1]); grid(1)
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1size_uharmonic7.pdf')


d(664458);figure(); 
d.dict['beam_area']=d.s1xgap*d.s1ygap;
subplot(2,2,1); d.plot('sum','fwhmx_mm','.g',hold=1); ylim([.1,.22]); grid(1)
subplot(2,2,2); d.plot('sum','fwhmy_mm','.g',hold=1); ylim([0.02,.1]); grid(1)
subplot(2,2,3); d.plot('beam_area','fwhmx_mm','.r',hold=1); ylim([.1,.22]); grid(1)
subplot(2,2,4); d.plot('beam_area','fwhmy_mm','.r',hold=1); ylim([0.02,.1]); grid(1)
subplots_adjust(wspace=.5, hspace=0.6)
savefig('/home/spc93/tmp/s1size_uharmonic11.pdf')


# focus vs time after opening FE shutter
figure();
d(664467)
subplot(2,2,1); d.plot('Elapsed','sum', hold=1); grid(1);
subplot(2,2,2); d.plot('Elapsed','T1dcmSi111', hold=1); grid(1);
subplot(2,2,3); d.plot('Elapsed','fwhmx_mm', hold=1); ylim([.15,.2]);grid(1);
subplot(2,2,4); d.plot('Elapsed','fwhmy_mm', hold=1); ylim([.06,.09]);grid(1);

savefig('/home/spc93/tmp/heatload.pdf')
figure()
d(664469)
subplot(2,2,1); d.plot('Elapsed','sum', hold=1); grid(1);
subplot(2,2,2); d.plot('Elapsed','T1dcmSi111', hold=1); grid(1);
subplot(2,2,3); d.plot('Elapsed','fwhmx_mm', hold=1); ylim([.16,.18]);grid(1);
subplot(2,2,4); d.plot('Elapsed','fwhmy_mm', hold=1); ylim([.06,.12]);grid(1);

savefig('/home/spc93/tmp/heatload9.pdf')
'''
