
#use Jython interpreter
#to dsplay: open h5 file/go to window for h5 file/ python/data/recip space/volume
#then slice data using isosurface

#https://confluence.diamond.ac.uk/display/I16/HKL+Mapping

import uk.ac.diamond.scisoft.analysis.diffraction.MillerSpaceMapper as Mapper
inpath='/dls/i16/data/2018/mt17842-1/'
outpath='/home/spc93/data/'



#scans=[670015, 670016]
#analyse separately
#for scan in scans:
#    Mapper.printCorners(inpath+'%06d.nxs' % scan, False)
#    #Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)

scan=678450
Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)


Mapper.processVolumeWithAutoBox('/dls/i16/data/2018/mt17842-1/678450_from_dat.nxs', '/home/spc93/data/678450.h5', "gaussian", 2.0, 2.0,True, 0.001)
