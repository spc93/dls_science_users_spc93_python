import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.max_rows',8)

path='/home/spc93/data/pdc_i16_sept_2017/'

#file='day14sept17_200_1.Spe'
#file='night_14sept17_200_1.Spe'
files=['night_14sept17_200_1.Spe', 'night_14sept17_200_2.Spe']
for file in files:
    print file
    #d=np.loadtxt(path+file).transpose()
    d=pd.read_csv((path+file),header=None)
    dat=d[0][12:-14]
    chans=array(range(4096))
    caldat=4118 #entrry for cal data
    c=float(str.split(str(d[0][caldat]))[0])
    m=float(str.split(str(d[0][caldat]))[1])
    en=m*chans+c

    plt.figure()
    plt.plot(en,dat)


'''


from dlstools import dataloader
path='/dls/i16/data/2016/cm14471-4/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')


samples=['diamond_6um','diamond_3um','silicon_nbs640']
for sample in samples:
    figure()
    for scan in range(613679, 613993+1):
        d(scan)
        if sample in d.note:
            if 'Bragg edge' in d.note:
                rat_edge=d.diode*1.0/d.IC2
            if 'no sample' in d.note:
                rat_null=d.diode*1.0/d.IC2
                double_ratio=rat_edge/rat_null
                subplot(2,2,1); plot(d.energy2, double_ratio, hold=1); axis('tight'); title('Sample = '+sample)
                subplot(2,2,2); plot(d.energy2-mean(d.energy2), double_ratio/mean(double_ratio), hold=1); axis('tight'); title('Sample = '+sample)
            if 'Bragg peak' in d.note:
                subplot(2,2,3); plot(d.energy2, d.APD, hold=1); axis('tight'); title('Sample = '+sample)
                subplot(2,2,4); plot(d.energy2-mean(d.energy2), d.APD*1.0/max(d.APD), hold=1); axis('tight'); title('Sample = '+sample)
    savefig('/home/spc93/tmp/%s.pdf' % sample)



#plot(d[0]*2,d[1]); title('TimePix coincidence spectrum');
#xlabel('Time (ns)'); ylabel('Counts')
#savefig('/home/spc93/tmp/bragedge1.pdf')

'''