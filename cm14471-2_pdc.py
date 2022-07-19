from dlstools import dataloader

path='/dls/i16/data/2016/cm14471-2/'
d=dataloader.dlsloader(path+'%i.dat')
close('all')



figure(); subplot(2,2,1)
d(591874)
hour=(d.TimeSec-d.TimeSec[0])/3600
d.plot(hour,'sum',hold=1); xlabel('hours')
d(591898)
subplot(2,2,2)
d.plot('mu','APD',hold=1);

savefig('/home/spc93/tmp/tmp.pdf')