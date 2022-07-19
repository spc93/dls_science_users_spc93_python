from dlstools import dataloader
close('all')

path='/dls/i16/data/2017/cm16772-1/'
d=dataloader.dlsloader(path+'%i.dat')

d(633777)
d.plot('idgap','ic1monitor');  
#ylim([1,max(d.ic1monitor)])
#gca().set_yscale('log')


d(633778)
d.plot('idgap','ic1monitor');  


d(633779)
d.plot('idgap','ic1monitor');  

#savefig('/home/i16user/tmp/tmp.pdf')
