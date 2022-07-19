from dlstools import dataloader


path='/dls/i16/data/2015/cm12169-3/'
d=dataloader.dlsloader(path+'%i.dat')
#p=dataloader.piloader(d, lambda obj: path+obj.pilatus100k_path_template)




print d(518482)
d.plot('kphi','sum')
d.plot('kphi','roiw_sum','r',hold=1);
d.plot('kphi',d.roiw_sum-d.sum,'r',hold=0);

#print p(2)
#figure(); pcolor(p.image0); axis('tight'); title(p.file)

#print d(447336)
#figure(); plot(d.eta, d.sum); title(d.file+'\n'+d.cmd); axis('tight'); grid(1)
#fit(pv_c)

