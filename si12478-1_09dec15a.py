from dlstools import dataloader
path='/dls/i10/data/2015/si12478-1/i10-'
d=dataloader.dlsloader(path+'%i.dat')

close('all')


subplot(2,2,1)
d(313771); plot(d.ips_field,1.0*d.macr19/d.macr16, hold=1); axis('tight')
d(313772); plot(d.ips_field,1.0*d.macr19/d.macr16,'r' ,hold=1); axis('tight'); grid(1); title('#'+str(d.datanumber)+' '+d.cmd)

subplot(2,2,2)
d(313783); plot(d.ips_field,1.0*d.macr19/d.macr16, hold=1); axis('tight')
d(313782); plot(d.ips_field,1.0*d.macr19/d.macr16,'r' ,hold=1); axis('tight'); grid(1); title('#'+str(d.datanumber)+' '+d.cmd)

