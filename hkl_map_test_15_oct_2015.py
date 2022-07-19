#module load dawn/nightly
print 'hello'

import uk.ac.diamond.scisoft.analysis.diffraction.MillerSpaceMapper as Mapper

#num, start hkl, vox size,

#Mapper.processVolume("/dls/i16/data/2015/cm12169-4/538039.nxs", "/scratch/temp/out2.h5", [50,32,24], [-0.10, -0.08, 8.94], 0.005)
#Mapper.processVolume("/dls/i16/data/2015/cm12169-4/538039.nxs", "/scratch/temp/out2_2.h5", [50*2,32*2,24*2], [-0.10, -0.08, 8.94], 0.0025)
#Mapper.processVolume("/dls/i16/data/2015/cm12169-4/538039.nxs", "/scratch/temp/out2_3.h5", [50*4,32*4,24*4], [-0.2, -0.2, 8.9], 0.0025)
#Mapper.processVolume("/dls/i16/data/2015/cm12169-4/538038.nxs", "/scratch/temp/out2_3.h5", [50,32,24], [-0.10, -0.08, 8.94], 0.005)