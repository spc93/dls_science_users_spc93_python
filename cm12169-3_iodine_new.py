from dlstools import dataloader


path='/dls/i16/data/2015/cm12169-3/'
d=dataloader.dlsloader(path+'%i.dat')
#p=dataloader.piloader(d, lambda obj: path+obj.pilatus100k_path_template)

## sample #1

d(525917)
iv=d.diode/d.ic1monitor
d(525919)
ih=d.diode/d.ic1monitor
d(525920)
i0=d.diode/d.ic1monitor
subplot(2,2,1); d.plot('energy2',-log(iv/i0), hold=1)
subplot(2,2,1); d.plot('energy2',-log(ih/i0), hold=1)
d(525924)
subplot(2,2,2); d.plot('chi','diode', hold=1)

d(525928)
iv=d.APD/d.ic1monitor
d(525930)
ih=d.APD/d.ic1monitor
d(525931)
i0=d.APD/d.ic1monitor
subplot(2,2,3); d.plot('energy2',-log(iv/i0), hold=1)
subplot(2,2,3); d.plot('energy2',-log(ih/i0), hold=1)
d(525933)
subplot(2,2,4); d.plot('chi','APD', hold=1)


#### sample 2
figure(2)
d(525944)
iv=d.APD/d.ic1monitor
d(525946)
ih=d.APD/d.ic1monitor
d(525947)
i0=d.APD/d.ic1monitor
subplot(2,2,1); d.plot('energy2',-log(iv/i0), hold=1)
subplot(2,2,1); d.plot('energy2',-log(ih/i0), hold=1)
d(525949)
subplot(2,2,2); d.plot('chi','APD', hold=1)


savefig('/home/i16user/tmp/tmp.pdf')


#d.plot('kphi','roiw_sum','r',hold=1);
#d.plot('kphi',d.roiw_sum-d.sum,'r',hold=0);

#print p(2)
#figure(); pcolor(p.image0); axis('tight'); title(p.file)

#print d(447336)
#figure(); plot(d.eta, d.sum); title(d.file+'\n'+d.cmd); axis('tight'); grid(1)
#fit(pv_c)

