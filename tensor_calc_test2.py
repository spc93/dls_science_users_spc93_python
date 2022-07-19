import TensorScatteringClass as ten
#from TensorScatteringClass import TensorScatteringClassMagrotExtension as ten

import numpy as np

t=ten.TensorScatteringClassMagrotExtension(CIFfile='/home/spc93/spc_cifs/CoCO3_icsd_61066.cif', Site='Co1');
#t.PlotIntensityInPolarizationChannelsVsMagrot('NonResMag', lam=1, hkl=np.array([0,0,9]), hkln=np.array([1,0,0]), K=2, Time=1, Parity=1, psideg=0, lk=np.array([.0,0,0]), sk=np.array([.1,0,0]), sigmapi='sigma')

#t.PlotIntensityInPolarizationChannelsVsMagrot('E1E1mag', lam=12.4/7, hkl=np.array([0,0,9]), hkln=np.array([1,0,0]), psideg=0, mk=np.array([1,0,0]), sigmapi='sigma')

#t.PlotIntensityVsPolarizationAnalyserRotation('NonResMag', lam = 12.4/5.22, hkl = np.array([0,0,9]), hkln = np.array([1,0,0]), psideg = 90, pol_eta_deg=range(361), K=None, Time=None, Parity=None, mk=None, lk=np.array([.0,0,0]), sk=np.array([.1,0,0]) , savefile=None)
#t.PlotIntensityVsPolarizationAnalyserRotation('Scalar', lam = 12.4/5.22, hkl = np.array([0,0,6]), hkln = np.array([1,0,0]), psideg = 90, pol_eta_deg=range(361), K=None, Time=None, Parity=None, mk=None, lk=np.array([.0,0,0]), sk=np.array([.1,0,0]) , savefile=None)



#t.PlotIntensityInPolarizationChannels('E2E2', lam=1, hkl=np.array([0,0,9]), hkln=np.array([1,0,0]), K=4, Time=1, Parity=1, sigmapi='sigma')


#(self, process, lam, hkl, hkln, psideg, pol_eta_deg, pol_th_deg = 45, stokesvec_swl = [0, 0, 1], K = None, Time = None, Parity = None, mk = None, lk = None, sk = None):

#process, lam, hkl, hkln, psideg, pol_eta_deg, pol_th_deg = 45, stokesvec_swl = [0, 0, 1], K = None, Time = None, Parity = None, mk = None, lk = None, sk = None, savefile=None):



t.PlotIntensityInPolarizationChannelsVsMagrot('NonResMag', lam=12.4/5.22, hkl=np.array([0,0,9]), hkln=np.array([1,0,0]), psideg=0, lk=np.array([0,0,0]), sk=np.array([1,0,0]), sigmapi='sigma')

