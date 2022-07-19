
#use Jython interpreter
#to dsplay: open h5 file/go to window for h5 file/ python/data/recip space/volume
#then slice data using isosurface

#https://confluence.diamond.ac.uk/display/I16/HKL+Mapping

import uk.ac.diamond.scisoft.analysis.diffraction.MillerSpaceMapper as Mapper
inpath='/dls/i16/data/2018/mt19259-1/'
#outpath='/home/spc93/data/'
outpath='/dls/science/users/spc93/processing/'

#scans=[670015, 670016]
#analyse separately
#for scan in scans:
#    Mapper.printCorners(inpath+'%06d.nxs' % scan, False)
#    #Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)

#scan=705081
#Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)
#Mapper.processVolumeWithAutoBox('/dls/i16/data/2018/mt17842-1/678450_from_dat.nxs', '/home/spc93/data/678450.h5', "gaussian", 2.0, 2.0,True, 0.001)


import time
#firstscan = 709533
firstscan = 709657
scan = firstscan
while True:
    try:
        print "=== Trying to map scan #%i" % scan
        Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)
        #Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.002)
        scan += 1
    except:
        print "=== Failed to map scan #%i - try again in five minutes" % scan
        time.sleep(300)