from dlstools import dataloader
from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/b16/data/2018/mt16638-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')

#scans=range(271665,271684+1,4) # before ajusting tomoPith
#scans=range(271710,271730+1,4) # after ajusting tomoPith
#scans=range(271730,271797+1,4) # after ajusting tomoPith
#scans=range(271799,271818+1,4) # after ajusting tomoPith
#scans=range(271819,271838+1,4) # after ajusting tomoPith
#scans=range(271839,271858+1,4) # after ajusting tomoPith
#scans=range(271959,271999+1,4) # after ajusting tomoPith
#scans=range(272307,272999+1,4) # after ajusting tomoPith


def PlotLoop(scans,fig):
    Iamscan=1
    figure(fig)
    while Iamscan==1:
        clf()
        subplot(2,2,1)
        for scan in scans:
            try:
                d(scan)
                plot(d.mlmXtal1Y,d.ai1,label='s1vo: '+str(d.s1_s1vo));
                xlabel('mlmXtal1Y');ylabel('ai1');
                print(scan)
            except:
                pass
            legend()
        subplot(2,2,2)
        for scan in scans:
            try:
                d(scan+1)
                plot(d.mlmXtal2Bragg,d.ai1,label='s1vo: '+str(d.s1_s1vo))
            except:
                pass
            xlabel('mlmXtal2Bragg');ylabel('ai1'); 
        subplot(2,2,3)
        for scan in scans:
            try:
                d(scan+2)
                plot(d.mlmXtal1Y,d.sum)
            except:
                pass
            xlabel('mlmXtal1Y');ylabel('Sum')
        
        subplot(2,2,4)
        for scan in scans:
            try:
                d(scan+3)
                plot(d.tomoTheta,d.fliproi_sum/d.ai1,label=str(scan+3))
            except:
                pass
            xlabel('TomoTheta');ylabel('ROI/ai1')
            legend()
        pause(5)
        return()
def PlotStatic(scans,fig):
    figure(fig)
    clf()
    subplot(2,2,1)
    for scan in scans:
        try:
            d(scan)
            plot(d.mlmXtal1Y,d.ai1,label='s1vo: '+str(d.s1_s1vo));
            xlabel('mlmXtal1Y');ylabel('ai1');
            print(scan)
        except:
            pass
        legend()
    subplot(2,2,2)
    for scan in scans:
        try:
            d(scan+1)
            plot(d.mlmXtal2Bragg,d.ai1,label='s1vo: '+str(d.s1_s1vo))
        except:
            pass
        xlabel('mlmXtal2Bragg');ylabel('ai1'); 
    subplot(2,2,3)
    for scan in scans:
        try:
            d(scan+2)
            plot(d.mlmXtal1Y,d.sum)
        except:
            pass
        xlabel('mlmXtal1Y');ylabel('Sum')
    
    subplot(2,2,4)
    for scan in scans:
        try:
            d(scan+3)
            plot(d.tomoTheta,d.fliproi_sum/d.ai1,label=str(scan+3))
        except:
            pass
        xlabel('TomoTheta');ylabel('ROI/ai1')
        legend()
    pause(5)
    return()

scans=range(272478,272507+1,4);
PlotStatic(scans,1)
scans=range(272551,272559+1,4);
PlotStatic(scans,1)
scans=range(272564,272570+1,4);
PlotLoop(scans,1)
show()
