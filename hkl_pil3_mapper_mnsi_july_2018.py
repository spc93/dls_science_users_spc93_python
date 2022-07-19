
#use Jython interpreter
#to dsplay: open h5 file/go to window for h5 file/ python/data/recip space/volume
#then slice data using isosurface

#https://confluence.diamond.ac.uk/display/I16/HKL+Mapping

#dataVis 100-1e9

import uk.ac.diamond.scisoft.analysis.diffraction.MillerSpaceMapper as Mapper
inpath='/dls/i16/data/2018/mt19259-1/'
outpath='/dls/science/users/spc93/processing/'

#27-30
scan=709532

print inpath+'%06d.nxs' % scan

Mapper.printCorners(inpath+'%06d.nxs' % scan, True)
Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)



