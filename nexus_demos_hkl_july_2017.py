import uk.ac.diamond.scisoft.analysis.diffraction.MillerSpaceMapper as Mapper

files=[ "/dls/i16/data/2017/mt16044-1/646%03d.nxs" % n for n in range(4,10)]
Mapper.processVolumeWithAutoBox(files, "/tmp/mnsi_fourpeaks.h5", "gaussian", 2.0, 2.0,True, 0.001)


#5e8
