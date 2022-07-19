from dlstools import dataloader
#from dlstools.quickfit import *
#from dlstools.dirty_fit import fit
path='/dls/i16/data/2017/mt14992-1/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')
C=pi/360/12.4 #multiply by delta*energy to get sin(th)/lam

##########
#various diffraction scans (P60 inlog book)
figure(); 
subplot(2,2,1); d(626889); d.plot('energy2','sum',hold=1); #SeK (bulk)
subplot(2,2,2); d(626891); d.plot('energy2','sum',hold=1); #BiL3 (bulk)
subplot(2,2,3); d(626897); d.plot('energy2','sum',hold=1); #BiL1 (bulk)

figure(); 
subplot(2,2,1); d(626915); d.plot('eta','sum',hold=1); #Al2O3 0,0,18
subplot(2,2,2); d(626928); d.plot('delta','roi2_sum',hold=1);gca().set_yscale('log'); #00l refs from film
subplot(2,2,3); d(626930); d.plot('kphi','roi2_sum',hold=1);gca().set_yscale('log'); #00l refs from film
subplot(2,2,4); d(626969); d.plot('alpha','roi2_sum',hold=1);gca().set_yscale('log'); #00l refs from film

###########
#tine scans with different ppp offsets - small offsets show more jitter (P58) 
figure()
iscans=range(627014,627032+1)
for iscan in iscans:
    d(iscan); d.plot('x','ic1', hold=1)
############
# diffraction scans, +ve offset P62

figure(); 

subplot(2,2,1); d(627070); d.plot('l',d.surf_sum/d.ic1monitor,hold=1); #
subplot(2,2,1); d(627091); d.plot('l',d.surf_sum/d.ic1monitor,hold=1); #

figure(); 
for irun in range(627107, 627113):
    d(irun); d.plot('l',d.fracdiff,hold=1);
    d.plot('l',d.norm_sum*1e-9,hold=1);
ylim([-.001,0.002])

#########################
#short scans range -ve offset P64
figure(); 
d(627127); sumrat=d.fracdiff*0; count=0
for irun in range(627127, 627134):
    d(irun); d.plot('l',d.fracdiff,hold=1);
    d.plot('l',d.norm_sum*1e-9,hold=1);
    sumrat+=d.fracdiff; count+=1;
ylim([-.001,0.002]); grid(1)
#mean value
figure(); plot(d.l, sumrat/count)
ylim([-.001,0.002]); grid(1)

############################
#reflectivity with diode
figure()
d(627145); d.plot(d.delta*d.en*C,'diode',hold=1); xlabel('sin($\\theta$)/$\lambda$')
ylim([0,10])
qq=(0.01026-0.005)/10
dd=1./2/qq  #film thickness (A) (950.6)

########################
#difference measurement (reflectivity) - ppp offset zero (by accident)
figure()
d(627155); d.plot('delta',d.fracdiff,hold=1)
d.plot('delta',d.norm_sum*1e-7,hold=1)
ylim([-.002,0.020])

################
# diffraction scans, +ve offset P67
figure(); 
d(627160); sumrat=d.fracdiff*0; count=0
for irun in range(627160, 627166+1):
    d(irun); d.plot('l',d.fracdiff,hold=1);
    d.plot('l',d.norm_sum*1e-9,hold=1);
    sumrat+=d.fracdiff; count+=1;
ylim([-.001,0.001]); grid(1)

#################
#ion chamber noise - analogue in and versus V-F (P68)
d(627169)
figure(); d.plot('x',d.ic1monitor/max(d.ic1monitor),hold=1)
d.plot('x',d.IC2*1.0/max(d.IC2),hold=1)

#################
#reflectivity: circ +ve, circ -v3, linear
figure()
d(627190); 
#d.plot('delta',d.fracdiff,'r',hold=1); #circpol offset +ve
d.plot(d.delta*d.en*C,d.fracdiff,'r',hold=1); #circpol offset +ve
#d.plot('delta',d.norm_sum*1e-7,hold=1)
d.plot(d.delta*d.en*C,d.norm_sum*1e-7,hold=1)
d(627191); 
#d.plot('delta',d.fracdiff,'b',hold=1) #circpol offset -ve
d.plot(d.delta*d.en*C,d.fracdiff,'b',hold=1) #circpol offset -ve
#d(627192); d.plot('delta',d.fracdiff,'g',hold=1) #linear after delta 0.65
ylim([-.01,0.040]); xlabel('sin($\\theta$)/$\lambda$'); ylabel('Fractional difference')

###################
#repeat diffraction scans but linear pol, const. Volts =4V (P70)
figure(); 
d(627205); sumrat=d.fracdiff*0; count=0
for irun in range(627205, 627210+1):
    d(irun); d.plot('l',d.fracdiff,hold=1);
    d.plot('l',d.norm_sum*1e-9,hold=1);
ylim([-.001,0.001]); grid(1)

#as above but no flipping - single pilatus measurement P70)
figure(); d(627211); d.plot('l','surf_sum'); 

#as above but using flipper with zero volts on digital out (no flipping)
figure(); #linear pol, const. Volts =4V, no volts on flipper
d(627212); sumrat=d.fracdiff*0; count=0
for irun in range(627212, 627213+1):
    d(irun); d.plot('l',d.fracdiff,hold=1);
    d.plot('l',d.norm_sum*1e-9,hold=1);
ylim([-.001,0.001]); grid(1)

###############################################
### new sample Bi2Te3 m2-0286
###############################################

#627234
#Te L3 edge fluo
##scan energy TeK-.02 TeK+.04 0.0005 pils 1
d(627234); figure(); d.plot('energy2','sum')
# 0,0,15 eta scan
d(627240);  figure(); d.plot('eta','roi2_sum') 

##day19feb17a.py
##627298 17:50 19/2/17
##pos x1 1
##for i in range(100):
##    pos y id(627240);  figure(); d.plot('eta','roi2_sum') 
##    for enval in [TeK, TeKhalf]:
##        pos energy enval
##        for hklval in [hh0015, hh006]:
##            pos hkl hklval
##            for ppoff in [-0.025, 0.025]:
##                pos x ppoff
##                pos ppb111 [enval, ppoff]
##                scancn hkl [0 0 .02] 41 checkbeam pil 1  mcroi
##                scancn hkl [0 0 .02] 41 checkbeam flipper12pilttl [1,1,8]
##TeK=4.341;TeKhalf=4.346

### 0 0 15 difference scans, on- and pre-edge, +ve and -ve circ pol (P74)
#some scans went to wrong hkl (??)
figure(); #0015
#ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]>10 and self.x<0 and "flipper12pilttl" in self.cmd')
ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]>10 and self.x<0 and "mcroi" in self.cmd')
subplot(2,2,1);
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', d.mcroi_sum*1e-7, hold=1); ylim(-.005, 0.01); 
ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]>10 and self.x<0 and "flipper12pilttl" in self.cmd')
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', 'fracdiff', hold=1); ylim(-.005, 0.01);  ylabel('fracdiff TeK -veoffset');

ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]>10 and self.x>0 and "mcroi" in self.cmd')
subplot(2,2,2);
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', d.mcroi_sum*1e-7, hold=1); ylim(-.005, 0.01);  
ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]>10 and self.x>0 and "flipper12pilttl" in self.cmd')
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', 'fracdiff', hold=1); ylim(-.005, 0.01);  ylabel('fracdiff TeK +veoffset');

ff=d.findscans(range(627298,627340+1),'self.en>4.343 and self.l[0]>10 and self.x<0 and "mcroi" in self.cmd')
subplot(2,2,3);
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', d.mcroi_sum*1e-7, hold=1); ylim(-.005, 0.01);  ylabel('fracdiff TeKhalf -veoffset');
ff=d.findscans(range(627298,627340+1),'self.en>4.343 and self.l[0]>10 and self.x<0 and "flipper12pilttl" in self.cmd')
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', 'fracdiff', hold=1); ylim(-.005, 0.01);   ylabel('fracdiff TeKhalf -veoffset');

ff=d.findscans(range(627298,627340+1),'self.en>4.343 and self.l[0]>10 and self.x>0 and "mcroi" in self.cmd')
subplot(2,2,4);
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', d.mcroi_sum*1e-7, hold=1); ylim(-.005, 0.01);  ylabel('fracdiff TeKhalf +veoffset');
ff=d.findscans(range(627298,627340+1),'self.en>4.343 and self.l[0]>10 and self.x>0 and "flipper12pilttl" in self.cmd')
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', 'fracdiff', hold=1); ylim(-.005, 0.01); ylabel('fracdiff TeKhalf +veoffset');

### 0 0 6 difference scans, on- and pre-edge, +ve and -ve circ pol (P75)
figure();
#ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]>10 and self.x<0 and "flipper12pilttl" in self.cmd')
ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]<10 and self.x<0 and "mcroi" in self.cmd')
subplot(2,2,1);
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', d.mcroi_sum*1e-7, hold=1); 
ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]<10 and self.x<0 and "flipper12pilttl" in self.cmd')
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', 'fracdiff', hold=1); ylim(-.005, 0.02);  ylabel('fracdiff TeK -veoffset');

ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]<10 and self.x>0 and "mcroi" in self.cmd')
subplot(2,2,2);
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', d.mcroi_sum*1e-7, hold=1);
ff=d.findscans(range(627298,627340+1),'self.en<4.343 and self.l[0]<10 and self.x>0 and "flipper12pilttl" in self.cmd')
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', 'fracdiff', hold=1); ylim(-.005, 0.02);  ylabel('fracdiff TeK +veoffset');

ff=d.findscans(range(627298,627340+1),'self.en>4.343 and self.l[0]<10 and self.x<0 and "mcroi" in self.cmd')
subplot(2,2,3);
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', d.mcroi_sum*1e-7, hold=1); ylabel('fracdiff TeKhalf -veoffset');
ff=d.findscans(range(627298,627340+1),'self.en>4.343 and self.l[0]<10 and self.x<0 and "flipper12pilttl" in self.cmd')
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', 'fracdiff', hold=1); ylim(-.005, 0.03);   ylabel('fracdiff TeKhalf -veoffset');

ff=d.findscans(range(627298,627340+1),'self.en>4.343 and self.l[0]<10 and self.x>0 and "mcroi" in self.cmd')
subplot(2,2,4);
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', d.mcroi_sum*1e-7, hold=1);  ylabel('fracdiff TeKhalf +veoffset');
ff=d.findscans(range(627298,627340+1),'self.en>4.343 and self.l[0]<10 and self.x>0 and "flipper12pilttl" in self.cmd')
for scan in ff:
    d(scan); d.plot('l',d.mcroi_sum, hold=1);
    d.plot('l', 'fracdiff', hold=1); ylim(-.005, 0.03); ylabel('fracdiff TeKhalf +veoffset');

###############################################
### new sample Sb2Te3 m2-057 (very old sample used for early I16 tests)
###############################################

#last scans - old Te sample (new ones showed no fringes!)
#day19feb17a.py
#627391-99
#pos x1 1
#for i in range(100):
#    pos y i
#    for enval in [TeK, TeKhalf]:
#        pos energy enval
#        for ppoff in [-0.025, 0.025]:
#            pos x ppoff
#            pos ppb111 [enval, ppoff]
#            scan delta 1 3 .005 eta .5+thoff .0025 t 1 flipper12apdttlvf [1,1,4]
#
#series of straight reflectivity scans, showing intensity increasing with time!
normalscans=range(627393,627398+1)
figure()
for scan in normalscans:
    d(scan)
    d.plot('delta','APD',hold=1,label=str(d.date))
legend()

#flipper scans
figure()
subplot(2,2,1)
d(627393); 
#d.plot('delta','fracdiff', hold=1)
d.plot(d.delta*d.en*C,'fracdiff', hold=1)
#d.plot('delta',d.APD*5e-8, hold=1)
d.plot(d.delta*d.en*C,d.APD*5e-8, hold=1)
d(627397); 
#d.plot('delta','fracdiff', hold=1); 
d.plot(d.delta*d.en*C,'fracdiff', hold=1);
title('Te pre-edge -ve pol.'); xlabel('sin($\\theta$)/$\lambda$')
ylim([-.02,0.1])

subplot(2,2,2)
d(627394); 
#d.plot('delta','fracdiff', hold=1)
d.plot(d.delta*d.en*C,'fracdiff', hold=1)
#d.plot('delta',d.APD*2e-8, hold=1)
d.plot(d.delta*d.en*C,d.APD*2e-8, hold=1)
d(627398); 
#d.plot('delta','fracdiff', hold=1); 
d.plot(d.delta*d.en*C,'fracdiff', hold=1); 
title('Te pre-edge +ve pol.'); xlabel('sin($\\theta$)/$\lambda$')
ylim([-.02,0.1])

subplot(2,2,3)
d(627395); 
#d.plot('delta',d.APD*2e-8, hold=1);  
d.plot(d.delta*d.en*C,d.APD*2e-8, hold=1);  
#d.plot('delta','fracdiff', hold=1)
d.plot(d.delta*d.en*C,'fracdiff', hold=1);  
title('Te main-edge -ve pol.'); xlabel('sin($\\theta$)/$\lambda$')
ylim([-.02,0.1])

subplot(2,2,4)
d(627396); 
#d.plot('delta',d.APD*2e-8, hold=1); 
d.plot(d.delta*d.en*C,d.APD*2e-8, hold=1); 
#d.plot('delta','fracdiff', hold=1)
d.plot(d.delta*d.en*C,'fracdiff', hold=1)
title('Te main-edge +ve pol.'); xlabel('sin($\\theta$)/$\lambda$')
ylim([-.02,0.1])

qq=(0.00663-0.00368)/3
dd=1./2/qq  #film thickness (A) (508.47)

####
savefig('/home/spc93/tmp/plot.pdf')
