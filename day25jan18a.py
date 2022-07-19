from dlstools import dataloader
path='/dls/i16/data/2018/cm19668-1/'
d=dataloader.dlsloader(path+'%i.dat')
p=dataloader.tiffloader(d, lambda obj: path+obj.pilatus2m_path_template)
close('all')

crash

#figure()
#d(670820) #fluo
#d.plot('p2m_tresh', d.sum/max(d.sum),'g', hold=1)
#d(670824) #elast
#d.plot('p2m_tresh', d.p2mroi1_sum/max(d.p2mroi1_sum),'r', hold=1); grid(1)

#low gain - much better
#figure()
#d(670825) #elast
#d.plot('p2m_tresh', d.p2mroi1_sum/max(d.p2mroi1_sum),'g', hold=1)
#d(670830) #elast
#d.plot('p2m_tresh', d.p2mroi1_sum/max(d.p2mroi1_sum),'r', hold=1); grid(1)



#subplot(2,2,1); d.plot('energy2', 'ic1monitor', hold=1)
#d.plot('energy2', 'diode', hold=1)


d(670844)
maxpix=100
p2msum=[]
for p2m_num in d.path:
    print p2m_num
    p(p2m_num)
    #p2msum+=[sum(p.image_01)]
    p2msum+=[sum(p.image_01[p.image_01<maxpix])]
figure; 
p2msum=array(p2msum)
p2mnormsum=p2msum/d.ic1monitor
plot(d.energy2, p2mnormsum,'g', hold=1)


d(670845)
maxpix=1000
p2msum=[]
for p2m_num in d.path:
    print p2m_num
    p(p2m_num)
    #p2msum+=[sum(p.image_01)]
    p2msum+=[sum(p.image_01[p.image_01<maxpix])]
figure; 
p2msum=array(p2msum)
p2mnormsum=p2msum/d.ic1monitor/10
plot(d.energy2, p2mnormsum,'r', hold=1); grid(1)



savefig('/home/spc93/tmp/tmp.pdf')