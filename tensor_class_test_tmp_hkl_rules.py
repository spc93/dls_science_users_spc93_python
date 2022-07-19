
#import TensorScatteringClass as ten
import TensorScatteringClass_tmp_test_hkl_rules as ten


#t=ten.TensorScatteringClass(CIFfile='/home/spc93/spc_cifs/nchem.2848-s6.cif', Site='Cu1');
#t.PlotIntensityInPolarizationChannels('E1E1', lam=1, hkl=np.array([0,0,1]), hkln=np.array([1,0,0]), K=2, Time=1, Parity=1, mk=None, nk=None, sk=None, sigmapi='sigma')
#t=ten.TensorScatteringClass(CIFfile='/home/spc93/spc_cifs/GaFeO3_icsd_55840.cif', Site='Fe3');
#t.PlotIntensityInPolarizationChannels('E1E1', lam=1, hkl=np.array([0,0,3]), hkln=np.array([1,0,0]), K=2, Time=1, Parity=1, mk=None, lk=None, sk=None, sigmapi='sigma')
#t=ten.TensorScatteringClass(CIFfile='/home/spc93/spc_cifs/mnsi_icsd_673224.cif', Site='Mn1');
#t.PlotIntensityInPolarizationChannels('NonResmag', lam=1.8, hkl=np.array([0,0,6]), hkln=np.array([1,0,0]), lk=np.array([0,0,0]), sk=np.array([1,0,0]), sigmapi='sigma')
#t.TensorCalc(hkl=np.array([4,0,0]), K=1, Parity=+1, Time=+1)
######################### check why zero symops at site (min=1)

cifpath='/home/spc93/spc_cifs/'
ciffile=cifpath+'CoCO3_icsd_61066.cif'
t=ten.TensorScatteringClassMagrotExtension(CIFfile=ciffile, Site='Co1');
t.PlotIntensityInPolarizationChannels('E1E1', lam=1, hkl=np.array([0,0,9]), hkln=np.array([1,0,0]), K=2, Time=1, Parity=1, sigmapi='sigma')

Ftot = 0
for v in t._vec_list:
    #Ftot = 
    print v