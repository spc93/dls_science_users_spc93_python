from dlstools.dataloader import pdnx
from dlstools import dataloader
from nexusformat.nexus import *
import pandas as pd
dataloader._scandata_field='/entry1/pil100k'


p='/dls/i16/data/2017/cm16772-2/%i.nxs'


n=pdnx(p % 642326)
print n.scan
n.plot(x='sgome',y='sum')
print n.nx.tree
n.nx.plot()

# what's missing: default attribute to specify default MXdata to plot, and possibly axis attributes
# for I16 the detector images are not in the NeXus!
# standard field for scan 'scan:NXdata'?

# For discussion:
#    complete dump of values with GDA names as NXcollection (no meaning outside of beamline)
#    standard classes (NXinstrument etc) with links to raw dump (these do have well-defined semantic but are cumbersome to use interactively)
#    develop detailed and complete application definitions (multiple), adequate to ensure correct processing by automatic pipeline.




