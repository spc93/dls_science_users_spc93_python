


import sys
sys.path.append('/dls_sw/i16/software/python')
from dlstools.pdnx import *
from matplotlib.pyplot import *


# use new modified NeXus file with NXmx in a second subentry
p = '/dls/science/users/spc93/misc_nexus_data/modified/tmp/%i.nxs'
d = '/home/spc93/tmp/nexus/%i.dat'
e = '/home/spc93/tmp/nexus/%i.xls'


#n = pdnx(p % 815893) # load NeXus file into pdnx (pandas/NeXus) object 
n = pdnx(p % 815893, entry = None, data = None ) # load NeXus file into pdnx (pandas/NeXus) object 


n.nx['/entry1/scan'].keys()

entry =  getNexusSubentryWithDefinition(n.nx, definition = 'NXclassic_scan')
for key in n.nx[entry].keys():
    if 'NXdata' in str(type(n.nx[entry][key])):
        entrydata = '%s/%s' % (entry, key)
        _use_classicscan = True


entrydata

#n.to_srs(d % 815893)         # export to classic SRS .dat file
#n.to_srs_plus(d % 815894)    # export to enhanced SRS .dat file (metadata as key-value pair assignments)
#n.to_excel(e % 815894)       # export to spreadsheet



####################################
# does it work if entry and data are specified?
# if yes then why does it not work if not specified?
# don't spend much time on this if it works with entry and data
