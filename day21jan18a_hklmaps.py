
#use Jython interpreter
#to dsplay: open h5 file/go to window for h5 file/ python/data/recip space/volume
#then slice data using isosurface

#https://confluence.diamond.ac.uk/display/I16/HKL+Mapping

import uk.ac.diamond.scisoft.analysis.diffraction.MillerSpaceMapper as Mapper
inpath='/dls/i16/data/2018/mt17676-1/'
outpath='/home/spc93/data/'

#tt=tictoc()
#pos atten 5
#pos hkl hh22m2_pil
#T=26.5, I=0V; #670015
#scancn hkl [.0001, .0001, -.0001] 401 pil 2
#print tt
#back to kepco. field on for conical peaks
#pos x22_anout -2.5*2.0
#pos hkl hh22m2_pil
#670016
#scancn hkl [.0001, .0001, -.0001] 401 pil 2

scans=[670015, 670016]
#analyse separately
for scan in scans:
    Mapper.printCorners(inpath+'%06d.nxs' % scan, False)
    #Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)




#B=0.2 T (permanent magnet)
#pos do do.pil
#pos delta 0
#pos atten 5
#670506
#for tval in frange(22,28,1):
#    pos tset tval
#    pos w 120
#    pos phi -18
#    pos hkl hh22m2_pil
#    scancn hkl [.0001, .0001, -.0001] 401 tset Ta Tb euler pil 2


scans=range(670506,670506+8)
#analyse separately
for scan in scans:
    Mapper.printCorners(inpath+'%06d.nxs' % scan, False)
    #Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)
