#con-alt-enter

#use Jython interpreter
#to dsplay: open h5 file/go to window for h5 file/ python/data/recip space/volume
#then slice data using isosurface

#https://confluence.diamond.ac.uk/display/I16/HKL+Mapping

import uk.ac.diamond.scisoft.analysis.diffraction.MillerSpaceMapper as Mapper
from time import sleep

outpath='/dls/science/users/spc93/processing/'

inpath='/dls/i16/data/2019/cm22970-1/'
scans=[736956]

for scan in scans:
#    #Mapper.printCorners(inpath+'%06d.nxs' % scan, False)
    Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)





















