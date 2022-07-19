from dlstools import dataloader
close('all')

path='/dls/i16/data/2017/nt17083-1/'
d=dataloader.dlsloader(path+'%i.dat')

figure(); 
d(633565); d.plot('energy2','roi1_sum', hold=1)
d(633566); d.plot('energy2','roi1_sum', hold=1)
d(633567); d.plot('energy2','roi1_sum', hold=1)


figure(); 
d(633672); d.plot('energy2','roi1_sum', hold=1)


'''
#overnight horizontal drift scans
scans=arange(499391,500205)
centres=[]; widths=[]; timesec=[]
for scan in scans:
    d(scan)
    figure(42);hold(0); plot(d.sx[1:], diff(d.diode)); axis('tight'); 
    fit(gauss_c);
    centres+=[gauss_c.p[1]]
    widths+=[gauss_c.p[2]]
    timesec+=[d.TimeSec]
centres=array(centres)
hdrift=(centres-centres[0])*1000.0
widths=array(widths)*1000.0
timesec=array(timesec); 
hours=(timesec-timesec[0])/3600.0


figure(10); 
subplot(2,2,1); plot(hours, hdrift); xlabel('Hours'); ylabel('Horiz. pos. (um)');  hold(0); grid('on'); axis('tight')
subplot(2,2,2); plot(hours, widths); xlabel('Hours'); ylabel('FWHM (um)'); hold(0); grid('on'); axis('tight')


#horizontal focus vs vertical position on mirror
scans=arange(500219,500255)
centres=[]; widths=[]; timesec=[]; vmtrans=[];
for scan in scans:
    d(scan)
    figure(42);hold(0); plot(d.sx[1:], diff(d.diode)); axis('tight'); 
    fit(gauss_c);
    centres+=[gauss_c.p[1]]
    widths+=[gauss_c.p[2]]
    timesec+=[d.TimeSec]
    vmtrans+=[d.kbm1_y]
centres=array(centres)
hdrift=(centres-centres[0])*1000.0
widths=array(widths)*1000.0
timesec=array(timesec); 
hours=(timesec-timesec[0])/3600.0
vmtrans=array(vmtrans)

figure(10); 
subplot(2,2,3); plot(vmtrans, hdrift); xlabel('vmtrans'); ylabel('Horiz. pos. (um)');  hold(0); grid('on'); axis('tight')
subplot(2,2,4); plot(vmtrans, widths); xlabel('vmtrans'); ylabel('FWHM (um)'); hold(0); grid('on'); axis('tight')

#horizontal position vs vtrans
scans=arange(500279,500296)
centres=[]; widths=[]; timesec=[]; ssx=[]
for scan in scans:
    d(scan)
    figure(42);hold(0); plot(d.sx[1:], diff(d.diode)); 
    plot(d.sx[23:45], diff(d.diode[22:45])); ######### zoomed in for better fit
    fit(gauss_c);
    centres+=[gauss_c.p[1]]
    widths+=[gauss_c.p[2]]
    timesec+=[d.TimeSec]
    ssx+=[d.s5xtrans]
centres=array(centres)
hdrift=(centres-centres[0])*1000.0
widths=array(widths)*1000.0
timesec=array(timesec); 
hours=(timesec-timesec[0])/3600.0
ssx=array(ssx)

figure(10); 
subplot(2,2,1); plot(ssx, hdrift); xlabel('horiz beam pos (mm)'); ylabel('Horiz. pos. (um)');  hold(0); grid('on'); axis('tight')
subplot(2,2,2); plot(ssx, widths); xlabel('horiz beam pos (mm)'); ylabel('FWHM (um)'); hold(0); grid('on'); axis('tight')


d(500331); #long overnight drift scan
figure(10); 
hours=(d.TimeSec-d.TimeSec[0])/3600.0
vdrift=(d.peakx_mm-d.peakx_mm[0])*1000
hdrift=(d.peaky_mm-d.peaky_mm[0])*1000
subplot(2,2,1); plot(hours, hdrift); xlabel('Hours'); ylabel('Horiz. pos. (um)');  hold(0); grid('on'); axis('tight'); ylim([-5,10]); title('#'+str(d.datanumber))
subplot(2,2,2); plot(hours, vdrift); xlabel('Hours'); ylabel('Vert. pos. (um)');  hold(0); grid('on'); axis('tight');  ylim([-15,0]);
subplot(2,2,3); plot(hours, d.fwhmy_mm*1000); xlabel('Hours'); ylabel('Horiz. fwhm. (um)');  hold(0); grid('on'); axis('tight'); ylim([0,50])
subplot(2,2,4); plot(hours, d.fwhmx_mm*1000); xlabel('Hours'); ylabel('Vert. fwhm. (um)');  hold(0); grid('on'); axis('tight'); ylim([0,50])
savefig('/home/i16user/tmp/tmp.pdf')

d(500333); #high-speed short scan
figure(10); 
hours=(d.TimeSec-d.TimeSec[0])/3600.0
vdrift=(d.peakx_mm-d.peakx_mm[0])*1000
hdrift=(d.peaky_mm-d.peaky_mm[0])*1000
subplot(2,2,1); plot(hours, hdrift); xlabel('Hours'); ylabel('Horiz. pos. (um)');  hold(0); grid('on'); axis('tight'); ylim([-0.5,1]); title('#'+str(d.datanumber))
subplot(2,2,2); plot(hours, vdrift); xlabel('Hours'); ylabel('Vert. pos. (um)');  hold(0); grid('on'); axis('tight'); ylim([-.5,1])
subplot(2,2,3); plot(hours, d.fwhmy_mm*1000); xlabel('Hours'); ylabel('Horiz. fwhm. (um)');  hold(0); grid('on'); axis('tight'); ylim([20,25])
subplot(2,2,4); plot(hours, d.fwhmx_mm*1000); xlabel('Hours'); ylabel('Vert. fwhm. (um)');  hold(0); grid('on'); axis('tight'); ylim([20,25])
savefig('/home/i16user/tmp/tmp.pdf')

d(500335); #high-speed short scan; move hmtrans in/out by 1 mm 
figure(10); 
hours=(d.TimeSec-d.TimeSec[0])/3600.0
vdrift=(d.peakx_mm-d.peakx_mm[0])*1000
hdrift=(d.peaky_mm-d.peaky_mm[0])*1000
subplot(2,2,1); plot(hours, hdrift); xlabel('Hours'); ylabel('Horiz. pos. (um)');  hold(0); grid('on'); axis('tight'); ylim([-0.5,1]); title('#'+str(d.datanumber))
subplot(2,2,2); plot(hours, vdrift); xlabel('Hours'); ylabel('Vert. pos. (um)');  hold(0); grid('on'); axis('tight'); ylim([-.5,1])
subplot(2,2,3); plot(hours, d.fwhmy_mm*1000); xlabel('Hours'); ylabel('Horiz. fwhm. (um)');  hold(0); grid('on'); axis('tight'); ylim([20,25])
subplot(2,2,4); plot(hours, d.fwhmx_mm*1000); xlabel('Hours'); ylabel('Vert. fwhm. (um)');  hold(0); grid('on'); axis('tight'); ylim([20,25])
savefig('/home/i16user/tmp/tmp.pdf')
'''

#savefig('/home/i16user/tmp/tmp.pdf')
