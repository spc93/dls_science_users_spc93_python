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
#scans=range(272150,272214+1,5) # after ajusting tomoPith
#scans=range(272478,272507+1,5) # after ajusting tomoPith
scans=range(272551,272607+1,5) # after ajusting tomoPith
Iamscan=1

while Iamscan==1:
    clf()
    subplot(2,3,1)
    for scan in scans:
        try:
            d(scan)
            plot(d.mlmXtal1Y,d.ai1,label='s1vo: '+str(d.s1_s1vo)+'_#'+str(scan));
            xlabel('mlmXtal1Y');ylabel('ai1');
            #print(scan)
        except:
            pass
        legend()
    subplot(2,3,2)
    for scan in scans:
        try:
            d(scan+1)
            plot(d.mlmXtal2Bragg,d.ai1,label='s1vo: '+str(d.s1_s1vo)+'_#'+str(scan+1))
        except:
            pass
        xlabel('mlmXtal2Bragg');ylabel('ai1'); 
        legend()
    subplot(2,3,3)
    for scan in scans:
        try:
            d(scan+2)
            plot(d.mlmXtal1Y,d.sum,label='_#'+str(scan+2))
        except:
            pass
        xlabel('mlmXtal1Y');ylabel('Sum')
        legend()
    subplot(2,3,4)
    for scan in scans:
        try:
            d(scan+3)
            plot(d.tomoTheta,d.fliproi_sum/d.ai1,label='_#'+str(scan+3))
        except:
            pass
        xlabel('TomoTheta');ylabel('ROI/ai1')
        legend()
    subplot(2,3,5)
    for scan in scans:
        try:
            d(scan+4)
            plot(d.x,d.fracdiff,label='s1vo:'+str(d.s1_s1vo)+'_#'+str(scan+4))
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
            s1xvo+=[mean(d.s1vo)]
            fracdiff+=[mean(d.fracdiff)]
            err+=[std(d.fracdiff)]
            errorbar(s1xvo, fracdiff,err,label='#'+str(scan+4));
            xlabel('Slit1 Height');ylabel('norm_diff');
        except:
            pass
        legend()
    grid(1)
    pause(10)
    show()
    pause(5)
    show()




#fig = plt.figure()
#ax = fig.gca(projection='3d')
#d.sequence(e0scans,'psi')
#ax.plot_wireframe(d.s.energy2, d.s.psi, d.s.APD, cstride=0)
#axis('tight'); ax.set_zlim3d(0,2000)
#xlabel('Energy(keV)'); ylabel('Azimuth');
#title('300 Energy scans pi-pi')
#savefig('/home/spc93/tmp/tmp_e0.pdf')