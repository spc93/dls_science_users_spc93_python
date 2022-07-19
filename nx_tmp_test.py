#import sys
#sys.path.append('/dls_sw/i16/software/python')
#from dlstools.pdnx import *
#from matplotlib.pyplot import *
#%matplotlib notebook 
#p='/dls/i16/data/2019/cm22970-1/%i.nxs'
#n=pdnx(p % 734847)

import sys
sys.path.append('/dls_sw/i16/software/python')
from dlstools.pdnx import *
from matplotlib.pyplot import *
#%matplotlib notebook 

#re4ad soleil file with pdnx - fails
#file = '/dls/science/users/spc93/misc_nexus_data/soleil/MoS2_sample3_hkl_scan_00370.nxs'
#n=pdnx(file)

#read with nexpy - fails
import nexusformat.nexus as nxx
_nx = nxx.nxload(file,'r')
#_nx.nxsetencoding('latin_1')
nxx.nxsetencoding('latin_1')

#read with h5py - ok
#import h5py
#f = h5py.File(file, 'r')
#for key in f['exptest_06239'].keys():
#    print key
#    try:
#        for keyy in f['exptest_06239'][key]:
#            print keyy
#    except:
#        pass
    
  
    
#nexpy gui: module load python/ana then nexpy. same error
#seem to have current version
#copy file from usb again then copy to dropbox and ask Ray why it wont open    