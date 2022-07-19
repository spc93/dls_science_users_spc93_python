from dlstools import dataloader
from nexusformat.nexus import *
import matplotlib.pyplot as plt
import pandas as pd
import h5py
close('all')

datpath='/dls/i16/data/2017/cm16772-2/%i.dat'
p='/dls/i16/data/2017/cm16772-2/%i.nxs'



d=dataloader.dlsloader(datpath)
d(640927)



n=nxload(p % 640927)
print n.tree
n.entry1.scan_command
figure();plt.plot(array(n.entry1.default.idgap+0), array(n.entry1.default.ic1monitor+0))


c=nxload('/home/spc93/data/chopper.nxs')
c.entry.data.data.plot()
c['/entry/instrument/detector/polar_angle'].plot()

#use h5py directly
nn=h5py.File(p % 640927)

specfile=nxload('/home/spc93/data/spec_example.nxs')

nxload('/home/spc93/data/spec_example.nxs').plot()


p='/home/spc93/data/%i.nxs'

pd.set_option('display.max_rows',16,'display.width',1000)
_alias_dict={'data':'/entry1/default', '_repr':'/entry1/scan_command'}
class pynx2(h5py.File):
    def __init__(self,  filestr):
        self.filestr=filestr
        h5py.File.__init__(self, filestr)
        try:
            _data=pd.DataFrame(dict(self[_alias_dict['data']]))
            setattr(self, 'data', _data)
        except:
            pass
        self.repstr=filestr
        try:
            self.repstr+='\n'+str(self[_alias_dict['_repr']][:])
        except:
            pass          

    def __repr__(self):
        return str(self.repstr)

z=pynx2(p % 559068)

###### not working ###########
pd.set_option('display.max_rows',16,'display.width',1000)
_alias_dict={'data':'/entry1/default', '_repr':'/entry1/scan_command'}
class pynx3(tree.load):
    def __init__(self,  filestr):
        self.filestr=filestr
        nxload.__init__(self, filestr)
        try:
            _data=pd.DataFrame(dict(self[_alias_dict['data']]))
            setattr(self, 'data', _data)
        except:
            pass
        self.repstr=filestr
        try:
            self.repstr+='\n'+str(self[_alias_dict['_repr']][:])
        except:
            pass          

    def __repr__(self):
        return str(self.repstr)




