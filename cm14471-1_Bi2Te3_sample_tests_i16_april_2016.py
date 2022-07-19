from dlstools import dataloader

path='/dls/i16/data/2016/cm14471-1/'
d=dataloader.dlsloader(path+'%i.dat')

close('all')

subplot(2,2,1) #Bi2Te3 film [3 0 l] rod, alpha=0.5, 0.25, 1.0
d(586993); d.plot('l','roi1_sum','b', hold=1); 
d(586995); d.plot('l','roi1_sum','g', hold=1); 
d(586997); d.plot('l','roi1_sum','r', hold=1); 
title('Bi2Te3 film [3 0 L] rod\nalpha=0.5(b), 0.25(g), 1.0(r)')
gca().set_yscale('log'); axis('tight')

subplot(2,2,2) #Bi2Te3 film [1 1 l] rod, alpha=0.5, 0.25, 1.0
d(586994); d.plot('l','roi1_sum','b', hold=1); 
d(586996); d.plot('l','roi1_sum','g', hold=1); 
d(586998); d.plot('l','roi1_sum','r', hold=1); 
title('Bi2Te3 film [1 1 L] rod\nalpha=0.5(b), 0.25(g), 1.0(r)')
gca().set_yscale('log'); axis('tight')

savefig('/home/spc93/tmp/Bi2Te3 film.pdf')

figure();
subplot(2,2,1) #Bi2Te3 substrate [2 4 l] rod, alpha=0.5, 0.25, 1.0
d(586999); d.plot('l','roi1_sum','b', hold=1); d(587000); d.plot('l','roi1_sum','b', hold=1);
d(587005); d.plot('l','roi1_sum','g', hold=1); d(587006); d.plot('l','roi1_sum','g', hold=1);
d(587011); d.plot('l','roi1_sum','r', hold=1); d(587012); d.plot('l','roi1_sum','r', hold=1);
title('Bi2Te3 substrate [2 4 L] rod\nalpha=0.5(b), 0.25(g), 1.0(r)')
gca().set_yscale('log'); axis('tight')

subplot(2,2,2) #Bi2Te3 substrate [3 0 l] rod, alpha=0.5, 0.25, 1.0
d(587001); d.plot('l','roi1_sum','b', hold=1); d(587002); d.plot('l','roi1_sum','b', hold=1);
d(587007); d.plot('l','roi1_sum','g', hold=1); d(587008); d.plot('l','roi1_sum','g', hold=1);
d(587013); d.plot('l','roi1_sum','r', hold=1); d(587014); d.plot('l','roi1_sum','r', hold=1);
title('Bi2Te3 substrate [3 0 L] rod\nalpha=0.5(b), 0.25(g), 1.0(r)')
gca().set_yscale('log'); axis('tight')

subplot(2,2,3) #Bi2Te3 substrate [0 3 l] rod, alpha=0.5, 0.25, 1.0
d(587003); d.plot('l','roi1_sum','b', hold=1); d(587004); d.plot('l','roi1_sum','b', hold=1);
d(587009); d.plot('l','roi1_sum','g', hold=1); d(587010); d.plot('l','roi1_sum','g', hold=1);
d(587015); d.plot('l','roi1_sum','r', hold=1); d(587016); d.plot('l','roi1_sum','r', hold=1);
title('Bi2Te3 substrate [0 3 L] rod\nalpha=0.5(b), 0.25(g), 1.0(r)')
gca().set_yscale('log'); axis('tight')

savefig('/home/spc93/tmp/Bi2Te3 substrate.pdf')

#For Jython console
#Mapper.processVolume(String input, String output, String splitter, double hwhm, double scale, int mShape, double mStart, double... mDelta)
#Mapper.processVolume("/dls/i16/data/2015/mt12821-1/561400.nxs", "/tmp/out_tmp.h5", [200,200,70], [-0.05,-0.2,3.7], 0.001)

#Mapper.processVolume("/dls/i16/data/2016/cm14471-1/586993.nxs", "/tmp/586993.h5", [20,20,400], [3-0.2,-0.2,1.0], 0.02)

