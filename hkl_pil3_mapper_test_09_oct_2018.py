#con-alt-enter

#use Jython interpreter
#to dsplay: open h5 file/go to window for h5 file/ python/data/recip space/volume
#then slice data using isosurface

#https://confluence.diamond.ac.uk/display/I16/HKL+Mapping

import uk.ac.diamond.scisoft.analysis.diffraction.MillerSpaceMapper as Mapper
from time import sleep

#inpath='/dls/i16/data/2018/mt19842-1/'
#scans=[718527]
outpath='/dls/science/users/spc93/processing/'

#inpath='/dls/i16/data/2018/cm19668-5/'
#scans=[726466]

#inpath='/dls/i16/data/2019/cm22970-1/'
#scans=[736956]

inpath='/dls/i16/data/2019/cm22970-1/'
scans=[736956]

# these fail
#inpath = '/dls/science/users/spc93/misc_nexus_data/modified/tmp/'
#inpath = '/dls/science/users/spc93/misc_nexus_data/modified/'
#scans = [815893] # modified NeXus file with NXmx in subentry

for scan in scans:
#    #Mapper.printCorners(inpath+'%06d.nxs' % scan, False)
    Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)

#firstscan = 721299
#firstscan = 721644

#scan = firstscan
#while True:
#    try:
#        print "=== Attempt to process file: ", scan
#        Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)
#        scan+=1
#    except:
#        print "===  Failed to process file: ", scan
#        sleep(10)

# print "=== Attempt to process file: ", scan
# Mapper.processVolumeWithAutoBox(inpath+'%06d.nxs' % scan, outpath+'%06d.h5' % scan, "gaussian", 2.0, 2.0,True, 0.001)



















