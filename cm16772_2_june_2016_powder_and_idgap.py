from nexusformat.nexus import *
import matplotlib.pyplot as plt
import pandas as pd
close('all')

p='/dls/i16/data/2017/cm16772-2/%i.nxs'

n=nxload(p % 642621)

plt.figure()
plt.plot(array(n.entry1.default.energy2+0), array(n.entry1.default.C4+0))
