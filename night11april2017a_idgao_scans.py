from dlstools import dataloader
import pandas as pd
from scipy import interpolate
import lmfit

pd.set_option('display.max_rows',10)

close('all')


path='/dls/i16/data/2017/cm16772-1/'
d=dataloader.dlsloader(path+'%i.dat')


#d(633777); d.plot('idgap','ic1monitor')
d(633778); d.plot('idgap','ic1monitor')
#d(633779); d.plot('idgap','ic1monitor')

path='/dls_sw/i16/software/python/userscripts/steve/'
fields=pd.read_csv(path+'U27_GBH.DAT',sep='\t', usecols=[0,1], names=['Gap', 'Bmax'], skiprows=[0])
print fields

B=interpolate.interp1d(fields.Gap, fields.Bmax, kind='linear')#use quadratic?
Gap=interpolate.interp1d(fields.Bmax, fields.Gap, kind='linear')

def undulator_peak_energy(gap, harmonic, Ebeam, period):
    #period cm, Ebeam GeV, B T, period mm
    K=0.0934*period*B(gap)
    return 9.50*Ebeam**2/(1+K**2/2)/period*harmonic

def undulator_gap(Epeak, harmonic, Ebeam, period):
    K=np.sqrt(2*(9.50*Ebeam**2*harmonic/Epeak/period -1))
    B=K/0.0934/period
    #return K
    return Gap(B)

gapoffset=0.5
undulator_peak_energy(6.254+gapoffset, 5, 3.0, 27.0)
undulator_gap(4.0, 3, 3.0, 27.0)+gapoffset

energy=14.0
for uharm in range(1,20):
    try:
        g=undulator_gap(energy, uharm, 3.05, 27.0)-gapoffset
        e=undulator_peak_energy(g+gapoffset, uharm, 3.05, 27.0)
        print "=== Energy: %.1f Harmonic %i: Gap=%.1f Recalculated energy=%.4f" % (energy, uharm, g, e)
    except:
        pass


gapvec=np.array([8.117, 5.490, 7.040, 9.539, 5.071, 5.801, 6.695, 7.768, 9.292, 11.712])
harmvec=np.array([3, 9, 7, 5, 17, 15, 13, 11, 9, 7])
envec=np.array([4.0, 8.0, 8.0, 8.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0])


#def gap_errors(envec,harmvec,gapvec, ebeam, gapoffset, period):
#    return undulator_gap(envec, harmvec, ebeam, period)-gapoffset-gapvec


def gap_errors(en_harm_gap, ebeam, gapoffset, period):
    (envec, harmvec, gapvec)=en_harm_gap
    return undulator_gap(envec, harmvec, ebeam, period)-gapoffset-gapvec



params=Parameters()
params.add('ebeam', min=2.9, max=3.1, vary=True)
params.add('gapoffset', min=0.0 ,max=1.0, vary=True)
params.add('period', value=27.0, vary=False)

#minimize(gap_errors, params, args=(envec,harmvec,gapvec), method='leastsq')

gap_errors((envec,harmvec,gapvec), 3.05, 0.5, 27.0)


x=array([1,2,3,4,5])
y=array([0,1,4,7,2])

def gaussdiff((x, y), amp, cen, wid):
    return amp * exp(-(x-cen)**2 /wid)

gaussdiff((x, y), 4., 3., 1.)
pg=Parameters()
pg.add('amp', vary=True)
pg.add('cen', vary=True)
pg.add('wid', vary=True)


#minimize(gaussdiff, pg, args=(x,y))


######################## use lmfit to minimize gap_errors (look at lmfit examples) #####################

'''
kind : str or int, optional
        Specifies the kind of interpolation as a string
        ('linear', 'nearest', 'zero', 'slinear', 'quadratic, 'cubic'
        where 'slinear', 'quadratic' and 'cubic' refer to a spline
        interpolation of first, second or third order) or as an integer
        specifying the order of the spline interpolator to use.
        Default is 'linear'.


#fit only odd harmonics


idgap peaks
4 keV (harmonic 4,3)
[6.254,8.117]
8 keV
(harmonic 9,8,7,6,5,4,3)
[5.490, 6.221, 7.040, 8.162, 9.539, 11.898, 16.498]
14 keV
(harmonic 17,16,...5)
[5.071, 5.404, 5.801, 6.199, 6.695, 7.188, 7.768, 8.489, 9.292, 10.389, 11.712, 13.936, 17.925] 
'''