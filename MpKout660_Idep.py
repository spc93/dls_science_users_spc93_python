from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/b16/data/2018/mt16638-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')


scans=range(272133,272149+1,1) # after ajusting tomoPith
#scans=range(272214,272250+1,1) # after ajusting tomoPith
Iamscan=1

while Iamscan==1:
    clf()
    subplot(1,2,1)
    for scan in scans:
        try:
            d(scan)
            plot(d.x,d.fracdiff,label='s1vo:-0.45_#'+str(scan))
        except:
            pass
        xlabel('X');ylabel('Diff');
        legend()
    subplot(1,2,2)
    s1xvo=[];fracdiff=[]; err=[]
    for scan in scans:
        try:
            d(scan)
            s1xvo=[];fracdiff=[]; err=[]
            s1xvo+=[mean(d.magvolts)/2.5]
            fracdiff+=[mean(d.fracdiff)]
            err+=[std(d.fracdiff)]
            errorbar(s1xvo, fracdiff,err,label='#'+str(scan));
            xlabel('Magnet Current (Amps)');ylabel('norm_diff');
        except:
            pass
    grid(1)
    pause(10)
    show()
