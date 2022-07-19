from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/b16/data/2018/mt16638-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')


#scans=range(272323,274000+1,5) # after ajusting tomoPith
#scans=range(272374,274000+1,5) # after ajusting tomoPith
#scans=range(272643,272686+1,5) # after ajusting tomoPith
scans=range(272746,272846+1,5)
scans=range(272771-4,272846+1,5)
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
    subplot(2,3,2)
    for scan in scans:
        try:
            d(scan+1)
            plot(d.mlmXtal2Bragg,d.ai1,label='s1vo: '+str(d.s1_s1vo))
        except:
            pass
        legend()
        xlabel('mlmXtal2Bragg');ylabel('ai1'); 
    subplot(2,3,3)
    for scan in scans:
        try:
            d(scan+2)
            plot(d.mlmXtal1Y,d.sum)
        except:
            pass
        xlabel('mlmXtal1Y');ylabel('sum')
    
    subplot(2,3,4)
    for scan in scans:
        try:
            d(scan+3)
            plot(d.tomoTheta,d.fliproi_sum/d.ai1)
        except:
            pass
        xlabel('TomoTheta');ylabel('ROI/ai1')
    subplot(2,3,5)
    for scan in scans:
        try:
            d(scan+4)
            plot(d.x,d.fracdiff,label='s1vo: '+str(d.s1_s1vo)+'#'+str(scan+4))
        except:
            pass
        xlabel('X');ylabel('Diff');
        legend()
    subplot(2,3,6)
    s1xvo=[];fracdiff=[]; err=[]
    for scan in scans:
        try:
            d(scan+4)
            s1xvo=[];fracdiff=[]; err=[]
            s1xvo+=[mean(d.s1_s1vo)]
            fracdiff+=[mean(d.fracdiff)]
            err+=[std(d.fracdiff)]
            errorbar(s1xvo, fracdiff,err,label='#'+str(scan+4));xlabel('s1vo');ylabel('norm_diff');
        except:
            pass
    #print('At s1vo  : ')
    #print((s1xvo[0]))
    #print(fracdiff[0]/err[0])
    grid(1)
    pause(10)
    show()
