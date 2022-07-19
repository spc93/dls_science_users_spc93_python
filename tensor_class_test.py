import TensorScatteringClass as ten

#t=ten.TensorScatteringClass(CIFfile='/home/spc93/spc_cifs/nchem.2848-s6.cif', Site='Cu1');
#t.PlotIntensityInPolarizationChannels('E1E1', lam=1, hkl=np.array([0,0,1]), hkln=np.array([1,0,0]), K=2, Time=1, Parity=1, mk=None, nk=None, sk=None, sigmapi='sigma')
#t.print_tensors()




#t=ten.TensorScatteringClass(CIFfile='/home/spc93/spc_cifs/GaFeO3_icsd_55840.cif', Site='Fe3');
#t.PlotIntensityInPolarizationChannels('E1E1', lam=1, hkl=np.array([0,0,3]), hkln=np.array([1,0,0]), K=2, Time=1, Parity=1, mk=None, lk=None, sk=None, sigmapi='sigma')
#t.print_tensors()


#t=ten.TensorScatteringClass(CIFfile='/home/spc93/spc_cifs/mnsi_icsd_673224.cif', Site='Mn1');
#t.PlotIntensityInPolarizationChannels('NonResmag', lam=1.8, hkl=np.array([0,0,6]), hkln=np.array([1,0,0]), lk=np.array([0,0,0]), sk=np.array([1,0,0]), sigmapi='sigma')
#t.TensorCalc(hkl=np.array([4,0,0]), K=1, Parity=+1, Time=+1)
#t.print_tensors()



######################### check why zero symops at site (min=1)
cifpath='/home/spc93/spc_cifs/'
ciffile=cifpath+'CoCO3_icsd_61066.cif'


#mc = Cr.Crystal()
#mc.load_cif(ciffile)
#refs=mc.reflection_list(7.11)

t=ten.TensorScatteringClassMagrotExtension(CIFfile=ciffile, Site='Co1');
#t.PlotIntensityInPolarizationChannelsVsMagrot('NonResmag', lam=12.4/7, hkl=np.array([0,0,9]), hkln=np.array([1,0,0]), psideg=0, lk=np.array([0,0,0]), sk=np.array([1,0,0]), sigmapi='sigma')
#t.PlotIntensityInPolarizationChannelsVsMagrot('E1E1mag', lam=12.4/7, hkl=np.array([0,0,9]), hkln=np.array([1,0,0]), psideg=0, mk=np.array([1,0,0]), sigmapi='sigma')
t.PlotIntensityInPolarizationChannels('E2E2', lam=1, hkl=np.array([0,0,9]), hkln=np.array([1,0,0]), K=4, Time=1, Parity=1, sigmapi='sigma')

#t.PlotIntensityInPolarizationChannels('NonResMag', 12.4/5.22, array([0,0,9]), array([1,0,0]), psideg=None, K=None, Time=None, Parity=None, mk=None, lk=array([.1,0,0]), sk=array([1,0,0]), sigmapi='sigma', savefile=None)


