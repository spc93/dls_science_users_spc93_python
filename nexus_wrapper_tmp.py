from dlstools.dataloader import pdnx
#from dlstools.dataloader import dataloader
from dlstools import dataloader
from nexusformat.nexus import *
import pandas as pd

p='/dls/i16/data/2017/cm16772-2/%i.nxs'
pdat='/dls/i16/data/2017/cm16772-2/%i.dat'
d=dataloader.dlsloader(pdat)

dataloader._scandata_field='/entry1/default'
n=pdnx(p % 642621)
for i in range(642621, 642642):print pdnx(p % i).scan
n.plot(x='energy2',y=['C1','C2','C3','C4'], title=n.scan)

d(642652)


dataloader._scandata_field='/entry1/default'
n651=pdnx(p % 642651)
print n651.scan
n651.plot()

dataloader._scandata_field='/entry1/bpm'
n652=pdnx(p % 642652)
print n652.scan
n652.plot()

dataloader._scandata_field='/entry1/pil100k'
n325=pdnx(p % 642325)
print n325.scan
n325.plot()







#nx_scan_dict=dict(nx[_scandata_field])
#for key in nx_scan_dict.keys():
#    print "=== doing %s" % key
#    nx_scan_dict[key]=nx_scan_dict[key].nxdata
#df=pd.DataFrame(nx_scan_dict) 

#n=pdnx(p % 642652)

