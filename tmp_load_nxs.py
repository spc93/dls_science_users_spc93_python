from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit

path='/dls/i16/data/2014/cm4968-5/'
d=dataloader.dlsloader(path+'%i.dat')
p=dataloader.tiffloader(d, lambda obj: path+obj.pilatus100k_path_template)
#n=dataloader.nxsloader2(path+'%i.nxs')

d(482940)
plot(d.Energy, d.sum,'g')

#n(482940)
#plot(n.Energy.Energy[:], n.pil100k.sum[:],'ro')

#import scisoftpy
#a=scisoftpy.io.load(path+'482940.nxs')
#plot(a.entry1.instrument.Energy.Energy[:], a.entry1.instrument.pil100k.sum[:])


from scisoftpy.dictutils import ListDict
class ListDictLoader(ListDict):
    def __init__(self):
        ListDict.__init__(self)
        self.myattribute=42
    def hello(self):
        print '<><><> Hello <><><>'
ltest=ListDictLoader()
ltest['thing1']='cat'
ltest['thing2']='hat'
ltest.other='boo' #same as dict


### this removes only order dictionary items - thisis what we need
for odkey in ltest.__odict__:
    del  ltest[odkey] 
ltest.new='ha!'
