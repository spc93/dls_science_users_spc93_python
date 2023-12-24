import sys
sys.path.append('/dls_sw/apps/scisoftpy/2.7')
sys.path.append('/dls_sw/i16/software/python')
from dlstools import *
from dlstools.pdnx import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt


#p='/dls/i16/data/2022/cm28156-15/%i.nxs'
#n = pdnx(p % 926175)
#n.nx.plot()
#n

# easy to get quick plot but painful to plot other columns

#p='/dls/i21/data/2022/cm31147-3/i21-%i.nxs'
p1 = '/dls/i21/data/2022/cm31147-3/BNOO6/i21-%i.nxs'
#n = pdnx(p % 240978)
n = pdnx(p1 % 247219)
#n = pdnx(p1 % 247219, entry = '/entry', data ='/diamond_scan')

#n.nx.plot()
#plot(np.array(n.nx.entry.instrument.energy.value), np.array(n.nx.entry.instrument.diff1_c.data)); title('diff_1 fluo detector')
