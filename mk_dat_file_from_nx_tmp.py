import sys
sys.path.append('/dls_sw/i16/software/python')
from dlstools.pdnx import *
from matplotlib.pyplot import *



inpath = '/dls/i16/data/2020/mm23580-2/'
outpath = '/dls/science/users/spc93/misc_nexus_data/modified/'
filenum = 815893

p='/dls/science/users/spc93/misc_nexus_data/modified/%i.nxs'

#n=pdnx(p % 815893)
#n=pdnx(p % 815893, entry = '/entry1/scan/', data = None)
n=pdnx(p % 815893, entry = None, data = None)
#entry = _entry, measurement = _measurement


positioner_tree_list  =  n.nx.entry1.scan.positioners.tree.split('\n')
assignments_list  =  [assig.lstrip() for assig in positioner_tree_list if '=' in assig and not '@'in assig]
try:
    scan_command_assignment = ["scan_command = '%s'" % str(n.nx.entry1.scan.scan_command)]
except:
    scan_command_assignment = []
assignments_list  =  ['<MetaDataAtStart>'] + scan_command_assignment + assignments_list + ['</MetaDataAtStart>']


d = '/dls/science/users/spc93/misc_nexus_data/modified/%i.dat'

#extra_metadata = []
extra_metadata = assignments_list
n.to_csv(d % 815893, sep = '\t', index = False)
with open(d % 815893, 'r+') as f:
    content = f.read()
    f.seek(0, 0)
#    for headerline in n.nx.entry1.scan.measurement.scan_header:
    for headerline in list(n.nx.entry1.scan.measurement.scan_header) + extra_metadata:
        f.write(headerline + '\n')
    f.write(' &END\n')
    f.write(content)

#e = '/dls/science/users/spc93/misc_nexus_data/modified/%i.xls'
#n.to_excel(e % 815893)



####### test these in jupyter hub on laptop for demo

d = '/dls/science/users/spc93/misc_nexus_data/modified/%i.dat'
e = '/dls/science/users/spc93/misc_nexus_data/modified/%i.xls'

n.to_srs(d % 815894)
n.to_srs_plus(d % 815894)
n.to_excel(e % 815893)

#################

def trunc(flt):
    return float('%g' % flt)




### demo .dat file with %g and srs header
### second version with key pair values
### add to_srs and to_srs_ext methods

#def find(self, keystring=''):
#   for key_sequence in self.findkeys(keystring):
#      obj = self.nx
#      for key in key_sequence:
#         obj = obj[key]
#         print('.nx' + self._list_to_dot_sep_string(key_sequence) + ' : \t', obj)


