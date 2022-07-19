#import __main__
from __main__ import gca, plot, axis, xlim
from lmfit import Model
import numpy as np


def peak(xdat,ydat,nbgpts=1):
    '''
    [centre, width, sum, height, area, m, c]=peak(x,y,nbgpts=1)
    Returns centre, width, sum, height, area after subtracting a sloping
    linear background from the fist and nbgpts data points and the background parameters (m, c)
    The width is derived from the standard deviation, asuming a gaussian shape
    xdat, ydat are 1d arrays or lists
    '''
    x=np.array(xdat); 
    y=np.array(ydat); 
    npts=len(x); 
    xspan=x[-1]-x[0];                    #copy data into arrays
    m=(np.mean(y[npts-nbgpts:npts])-np.mean(y[0:nbgpts]))/(np.mean(x[npts-nbgpts:npts])-np.mean(x[0:nbgpts]))    #slope for linear b/g
    c=np.mean(y[0:nbgpts])-np.mean(x[0:nbgpts])*m;                            #intercept
    y=y-m*x-c;                                            #subtract background
    sumdat=sum(y);                                            #peak sum
    area=sumdat*xspan/npts;                                        #peak integral
    centre=sum(x*y)/sumdat;                                    #centroid calc.
    height=max(y);                                            #max y value after linear b/g
    fwhm_sd=np.sqrt(sum((x-centre)**2*y)/sumdat) * np.sqrt(8*np.log(2));  #fwhm from standard deviation (correct for gaussian)
    fwhm_area = area/height * 0.3989 *np.sqrt(8*np.log(2)) #fwhm from area/height (correct for gaussian) - more robust than sd with noisy data
    return [centre, fwhm_sd, fwhm_area, sumdat, height, area, m, c]    
 

### old fit functions   
#gauss_c=fit_func('gaussian + const',['area','centre','width','constant'],
#                    'area/width*np.sqrt(4*np.log(2)/np.pi)*np.exp(-4*np.log(2)*((x-centre)/width)**2)+constant',
#                    '[centre, width, sum, height, area, m, c, constant]=peak(x,y)')
#                    
#lor_c=fit_func('Lorentzian + const',['area','centre','width','constant'],                   
#                'area/width/(np.pi/2)/(1+4*((x-centre)/width)**2)+constant',
#                '[centre, width, sum, height, area, m, c, constant]=peak(x,y)')            
#                
#pv_c=fit_func('Pseudo-Voigt + const',['area','centre','width','lfrac','constant'],
#              'area/width/(lfrac*np.pi/2+(1-lfrac)*np.sqrt(np.pi/4/np.log(2)))*(lfrac/(1+4*((x-centre)/width)**2)+(1-lfrac)*np.exp(-4*np.log(2)*((x-centre)/width)**2))+constant',
#              '[centre, width, sum, height, area, m, c, constant]=peak(x,y); lfrac=1')  
####################



def gauss(x, area, cen, fwhm):
    return area/fwhm*np.sqrt(4*np.log(2)/np.pi)*np.exp(-4*np.log(2)*((x-cen)/fwhm)**2)

def poly2(x, m, c):
   return m * x + c

g_lin = Model(gauss) + Model(poly2)


#def gaussian(x, amp, cen, wid):
#   return amp * np.exp(-(x-cen)**2 / wid)
#gauss = Model(gaussian)
#
#def poly2(x, m, c):
#   return m * x + c
#lin_bg = Model(poly2)
#
#g_lin = gauss + lin_bg



def fit(func, aXis=None):
    '''
    fit(func): fit first line in current axis using lmfit Model instance func
    fit(func, axis): do the same for the line in a specified axis
    fit is designed to be simple with limited functionality
    Use lmfit directly for a more sophistcated fit
    Need to import pyplot first: from matplotlib.pyplot import *
    Need interactive graphics in Jupyter
    '''
    if aXis == None:
        aXis = gca()

    [xData, yData] = aXis.get_lines()[0].get_xydata().transpose()
    xL, xU = aXis.get_xlim()
    iI = (xData >= xL) & (xData <= xU)
    xData, yData = xData[iI], yData[iI]
    
    print 'xData', xData
    
    pk_prms = {}
    [pk_prms['cen'], fwhm_sd, pk_prms['fwhm'], pk_prms['sum'], pk_prms['amp'], pk_prms['area'], pk_prms['m'], pk_prms['c']] = peak(xData, yData)
    
    print 'params from peak:',pk_prms
    _prms = func.make_params()
    
    # assign fit parameter value to parameters that match the parameters from peak with others defaulting to zero
    for key in _prms.keys():
        try:
            _prms[key].value = pk_prms[key]
        except:
            _prms[key].value = 0
            
    #print xData
    #print yData
    #print pk_prms['wid']
    #print _prms['wid']
    result =  func.fit(yData, x=xData, params = _prms)
    #print result.params['wid']
    

    #for pname in result.params:
        #prm = result.params[pname]
	#print prm.name, prm.value, prm.stderr


    
    print result.fit_report()
 
    outstr = func.name+'\n'
    for pname in result.params:
        prm = result.params[pname]
	#_n_dec_places = int(-np.round(np.log10(self.params[name].stderr))) + 1
        try:
            _n_dec_places = max(0, int(-np.round(np.log10(prm.stderr))) + 1)
        except:
            _n_dec_places = 4 # in case parameter not varied then stderr is zero
	_fmt='%.' + str(_n_dec_places) + 'f'

        #print _fmt
        #self.out+='%-10s:  ' % name+str(self.params[name].value)+' +/- '+str(self.params[name].stderr)+'\n'
	outstr+='%10s:  ' % prm.name + '%15s' % (_fmt % prm.value)+' +/- '+ '%-10s' % str(_fmt % prm.stderr) +'\n'

    print outstr
 

    plot(xData, func.eval(x=xData, params=result.params),'r.'); axis('tight'); xlim(xL, xU); legend()
    #return result
    return outstr
    
# plot test data
#x = frange(20,30,.1)
#y = np.exp(-(x-25)**2) + rand(len(x))/10
#figure(); plot(x, y)
#fit(gauss)

