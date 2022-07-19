from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/b16/data/2018/mt16638-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')


#scans=range(272492,272600+1,1)+range(272585,272600+1,1) # after ajusting tomoPith
#scans=[272816,272822,272827,272832,272837,272844]
#scans=range(272849,272859+1,1);
scans=[272657,272662,272667,272672,272682,272863,272868,272873,272886,272896]

Iamscan=1

while Iamscan==1:
    clf()
    subplot(1,2,1)
    for scan in scans:
        try:
            d(scan)
            plot(d.x,d.fracdiff,label='-'+str(scan))
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
            s1xvo+=[mean(d.s1vo)]
            #s1xvo+=[mean(d.magvolts)/2.5]
            fracdiff+=[mean(d.fracdiff)]
            err+=[std(d.fracdiff)]
            errorbar(s1xvo, fracdiff,err,label='#'+str(scan));
            #xlabel('Magnetic Current(Amp)');ylabel('norm_diff');
            xlabel('s1vo');ylabel('norm_diff');
        except:
            pass
    grid(1)
    pause(10)
    show()
