##### get all Delta E values
##### re-measure CoCO3

from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit

close('all')

#outpath='/media/sf_Dropbox/tmp/'
outpath='/home/i16user/tmp/'

#expt='FeBO3'
#expt='MnCO3'
#expt='CoCO3c'
expt='NiCO3a' #first data sets
#expt='NiCO3b' #second data sets
#del dataset #comment out to avoid re-loading data

try:
    dataset
    assert expt==name
    print '=== dataset dict already created - skip'
except:
    
    if expt=='FeBO3':
        dataset=[{'pol':'sp', 'psi':24.34, 'res':'+', 'ymax':-1},
            {'pol':'ss', 'psi':24.34,'res':'+', 'ymax':-1},
            {'pol':'sp', 'psi':33.17,'res':'-','ymax':-1},
            {'pol':'ss', 'psi':33.17,'res':'-','ymax':30000},
            {'pol':'sp', 'psi':11,'res':'++','ymax':100000},
            {'pol':'ss', 'psi':11,'res':'++','ymax':20000},
            {'pol':'sp', 'psi':43,'res':'--','ymax':100000},
            {'pol':'ss', 'psi':43,'res':'--','ymax':100000},
            ]
        d=dataloader.dlsloader('/media/sf_data/FeBO3_I16_July_2014/%i.dat')
        first=455580
        last=456003
        delta=7.1125 #centre of E2 resonance
        
        for data in dataset: #get scan numbers
            data['scans']=d.findscans(range(first, last+1),'self.note=="magrot_scan_%s" and self.psi==%.2f and "magrot" in self.cmd' % (data['pol'],data['psi']))

    elif expt=='MnCO3':
        dataset=[{'pol':'sp', 'psi':20.0, 'res':'+', 'ymax':-1},
                 {'pol':'ss', 'psi':20.0,'res':'+', 'ymax':-1},
                 {'pol':'sp', 'psi':40.0,'res':'-','ymax':-1},
                 {'pol':'ss', 'psi':40.0,'res':'-','ymax':-1},
                 {'pol':'sp', 'psi':6.0,'res':'++','ymax':-1},
                 {'pol':'ss', 'psi':6.0,'res':'++','ymax':-1},
                 {'pol':'sp', 'psi':54.0,'res':'--','ymax':-1},
                 {'pol':'ss', 'psi':54.0,'res':'--','ymax':-1},
                 ]
        d=dataloader.dlsloader('/media/sf_data/mt9507-1/%i.dat')
        first=448809
        last=456003
        delta=6.5382 #centre of E2 resonance (from fluo)
        for data in dataset: #get scan numbers
            if data['pol']=='ss':
                pol='sigsig'
            elif data['pol']=='sp':
                pol='sigpi'
            else:
                pol=''
            data['scans']=d.findscans(range(first, last+1),'self.note=="magrot_scan_%s" and np.abs(self.psi-%.2f)< 1 and "magrot" in self.cmd' % (pol,data['psi']))

    elif expt=='CoCO3a':
        dataset=[{'pol':'sp', 'psi':50.0, 'res':'+', 'ymax':50000},
                 {'pol':'ss', 'psi':50.0,'res':'+', 'ymax':10000},
                 ]
        d=dataloader.dlsloader('/media/sf_data/mt8976-1/%i.dat')
        first=432454
        last=432636
        delta=7.70825 #centre of E2 resonance
        for data in dataset: #get scan numbers
            if data['pol']=='ss':
                polang=3
            elif data['pol']=='sp':
                polang=93
            else:
                pol=''
            data['scans']=d.findscans(range(first, last+1),'np.abs(self.stoke-%.2f)< 1 and np.abs(self.psi-%.2f)< 1 and "xps3motor1" in self.cmd' % (polang,data['psi']))

    elif expt=='CoCO3b':        #I16 Oct 2014 T=4K - large anisotropy
        dataset=[{'pol':'sp', 'psi':50.0, 'res':'--', 'ymax':-1},
                 {'pol':'ss', 'psi':50.0,'res':'--', 'ymax':-1},
                 {'pol':'sp', 'psi':33.4, 'res':'-', 'ymax':-1},
                 {'pol':'ss', 'psi':33.4,'res':'-', 'ymax':-1},
                 {'pol':'sp', 'psi':25.8, 'res':'+', 'ymax':-1},
                 {'pol':'ss', 'psi':25.8,'res':'+', 'ymax':-1},
                 ]
        d=dataloader.dlsloader('/media/sf_data/cm4968-5/%i.dat')
        first=482986;  last=483298
        delta=7.70850 #centre of E2 resonance
        for data in dataset: #get scan numbers
            data['scans']=d.findscans(range(first, last+1),'self.note=="magrot_scan_%s" and np.abs(self.psi-%.2f)< 1 and "magrot" in self.cmd' % (data['pol'],data['psi']))

    elif expt=='CoCO3c':        #I16 Oct 2014 T=13K - small anisotropy
        dataset=[{'pol':'sp', 'psi':10.7, 'res':'++', 'ymax':-1},
                 {'pol':'ss', 'psi':10.7,'res':'++', 'ymax':-1},
                 {'pol':'sp', 'psi':50.0, 'res':'--', 'ymax':-1},
                 {'pol':'ss', 'psi':50.0,'res':'--', 'ymax':-1},
                 {'pol':'sp', 'psi':33.4, 'res':'-', 'ymax':-1},
                 {'pol':'ss', 'psi':33.4,'res':'-', 'ymax':-1},
                 {'pol':'sp', 'psi':25.8, 'res':'+', 'ymax':-1},
                 {'pol':'ss', 'psi':25.8,'res':'+', 'ymax':-1},
                 ]
        d=dataloader.dlsloader('/media/sf_data/cm4968-5/%i.dat')
        first=483337;  last=483791
        delta=7.70850 #centre of E2 resonance
        for data in dataset: #get scan numbers
            data['scans']=d.findscans(range(first, last+1),'self.note=="magrot_scan_%s" and np.abs(self.psi-%.2f)< 1 and "magrot" in self.cmd' % (data['pol'],data['psi']))

    elif expt=='NiCO3a':        #I16 Dec 2014 T=5K 
        dataset=[{'pol':'sp', 'psi':7.8, 'res':'++', 'ymax':-1},
                 {'pol':'ss', 'psi':7.8,'res':'++', 'ymax':-1},
                 {'pol':'sp', 'psi':22.2, 'res':'+', 'ymax':-1},
                 {'pol':'ss', 'psi':22.2,'res':'+', 'ymax':-1},
                 {'pol':'sp', 'psi':28.4, 'res':'+', 'ymax':-1},
                 {'pol':'ss', 'psi':28.4,'res':'+', 'ymax':-1},
                 {'pol':'sp', 'psi':36.0, 'res':'-', 'ymax':-1},
                 {'pol':'ss', 'psi':36.0,'res':'-', 'ymax':-1},
                 ]
        d=dataloader.dlsloader('/dls/i16/data/2014/cm4968-5/%i.dat')
        first=488968;  last=489484
        delta=8.332 #centre of E2 resonance
        for data in dataset: #get scan numbers
            data['scans']=d.findscans(range(first, last+1),'self.note=="magrot_scan_%s" and np.abs(self.psi-%.2f)< 1 and "magrot" in self.cmd' % (data['pol'],data['psi']))
            
    elif expt=='NiCO3b':        #I16 Dec 2014 T=5K 
        dataset=[{'pol':'sp', 'psi':36.0, 'res':'-', 'ymax':-1},
                 {'pol':'ss', 'psi':36.0,'res':'-', 'ymax':-1},
                 {'pol':'sp', 'psi':40.9, 'res':'-', 'ymax':-1},
                 {'pol':'ss', 'psi':40.9,'res':'-', 'ymax':-1},
                 {'pol':'sp', 'psi':52.2, 'res':'-', 'ymax':-1},
                 {'pol':'ss', 'psi':52.2,'res':'-', 'ymax':-1},
                 {'pol':'sp', 'psi':68.1, 'res':'--', 'ymax':-1},
                 {'pol':'ss', 'psi':68.1,'res':'--', 'ymax':-1},
                 ]
        d=dataloader.dlsloader('/dls/i16/data/2014/cm4968-5/%i.dat')
        first=488968;  last=489756
        delta=8.332 #centre of E2 resonance
        for data in dataset: #get scan numbers
            data['scans']=d.findscans(range(first, last+1),'self.note=="magrot_scan_%s" and np.abs(self.psi-%.2f)< 1 and "magrot" in self.cmd' % (data['pol'],data['psi']))

name=expt
#plot_no=0
#nn=figure().number; #figure for polar plot

#for data in [dat for dat in dataset if dat['pol']=='sp']: #sp only
for ii in range(1,4+1):   figure(ii, figsize=(12, 12)).nextsubplot=1    #keep track of subplots with new figure attribute
for dat in dataset:
    if dat['pol']=='sp':
        polfig=1; fitfig=3;
    elif dat['pol']=='ss':
        polfig=2; fitfig=4;
    else:
        print '=== unkown pol'
        raise ValueError
    
    en=[]; fitp=[]; 
    #plot_no+=1
    for scan in dat['scans']:
        d(scan);
        
        if expt=='MnCO3' and dat['pol']=='ss':  #usually APD but two MPIX rois for MnCO3
            cts=d.mlroi_sum
        elif expt=='MnCO3' and dat['pol']=='sp':
            cts=d.rl90_sum
        else:
            cts=d.APD
        if expt=='CoCO3a':
            d.magrot=d.psi-d.xps3motor1 ### some incertainty over this...
            
        en+=[d.Energy];
        cos_sin_cos_sin_c.fit(d.magrot*pi/180, cts); fitp+=[cos_sin_cos_sin_c.p];
        figure(99); plot(d.magrot*pi/180, cts, d.magrot*pi/180,cos_sin_cos_sin_c(d.magrot*pi/180) ); axis('tight')
        #pause(0.1); # show plots updating
        
        if d.Energy<delta:
            colour='g'
        elif d.Energy==delta:
            colour='k'
        else:
            colour='r'
        cts[-1]=cts[0] #use 0 data instead of 360 deg to makes lies join
        titlestr='pol = %s psi = %.2f (%s)' % (dat['pol'], dat['psi'], dat['res'])
        figure(polfig); subplot(2,2, figure(polfig).nextsubplot, polar=True); plot(d.magrot*pi/180, cts*1./cts.max(), color=colour, linewidth=1.5 ); title(titlestr); suptitle(expt); 
        #if data['ymax']>0:  ylim(0,data['ymax'] )
    figure(polfig).nextsubplot+=1
    savefig('%s%s_%s_%s' % (outpath,expt,dat['pol'],'polar.pdf'))
    fitp=array(fitp);
    cosx_amp=fitp.transpose()[0]; sinx_amp=fitp.transpose()[1]; cos2x_amp=fitp.transpose()[2]; sin2x_amp=fitp.transpose()[3]; const=sin2x_amp=fitp.transpose()[4];
    figure(fitfig); subplot(2,2,figure(fitfig).nextsubplot); plot(en, cosx_amp, label='cosx_amp'); plot(en, sinx_amp, label='sinx_amp'); plot(en, cos2x_amp, label='cos2x_amp'); plot(en, sin2x_amp,label='sin2x_amp');  plot(en, const, label='const'); title('pol = '+dat['pol']+' psi = %.2f' % dat['psi']); grid(True) 
    figure(fitfig).nextsubplot+=1
    axis('tight'); title(titlestr);
    if dat['ymax']>0:  ylim(-dat['ymax']/3,dat['ymax'] )
    legend(loc='best'); suptitle(expt)
    savefig('%s%s_%s_%s' % (outpath,expt,dat['pol'],'fits.pdf'))
    
    print '=== pol=%s polfig=%i  fitfig=%i' % (dat['pol'], polfig, fitfig)
    for ii in range(1,4+1):   print '=== fig %i next subplot %i' % (ii, figure(ii).nextsubplot)


