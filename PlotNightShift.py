from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/b16/data/2018/mt16638-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')


scans=range(271879,272000+1,7) # after ajusting tomoPith
Iamscan=1

while Iamscan==1:
    clf()
    subplot(2,3,1)
    for scan in scans:
        try:
            d(scan)
            plot(d.mlmXtal1Y,d.ai1,label='s1vo: '+str(d.s1_s1vo));xlabel('mlmXtal1Y');ylabel('ai1');
            print(scan)
        except:
            pass
        legend()
    subplot(2,3,2)
    for scan in scans:
        try:
            d(scan+1)
            plot(d.mlmXtal2Bragg,d.ai1,label='s1vo: '+str(d.s1_s1vo))
            d(scan+3)
            plot(d.mlmXtal2Bragg,d.ai1,label='s1vo: '+str(d.s1_s1vo),ls='--')
        except:
            pass
        xlabel('mlmXtal2Bragg');ylabel('ai1'); 
    subplot(2,3,3)
    for scan in scans:
        try:
            d(scan+2)
            plot(d.mlmXtal1Y,d.sum)
            d(scan+4)
            plot(d.mlmXtal1Y,d.sum,ls='--')
        except:
            pass
        xlabel('mlmXtal1Y');ylabel('sum')
    
    subplot(2,3,4)
    for scan in scans:
        try:
            d(scan+5)
            plot(d.tomoTheta,d.fliproi_sum/d.ai1)
        except:
            pass
        xlabel('TomoTheta');ylabel('ROI/ai1')
    subplot(2,3,5)
    for scan in scans:
        try:
            d(scan+6)
            plot(d.x,d.fracdiff,label='s1vo: '+str(d.s1_s1vo)+'#'+str(scan+6))
        except:
            pass
        xlabel('X');ylabel('Diff');
        legend()
    subplot(2,3,6)
    flipruns=range(271885,271945,7)
    s1xvo=[];fracdiff=[]; err=[]
    for run in flipruns:
        print d(run)
        s1xvo+=[mean(d.s1_s1vo)]
        fracdiff+=[mean(d.fracdiff)]
        err+=[std(d.fracdiff)]
    errorbar(s1xvo, fracdiff,err,label='#'+str(flipruns));xlabel('s1vo');ylabel('fracdiff');
    grid(1)
    pause(10)
    show()
