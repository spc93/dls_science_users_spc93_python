import nexusformat.nexus as nx
import pandas as pd
import matplotlib
import numpy as np
from astropy.modeling.tests import data

pd.set_option('display.max_rows',8)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width',999)

#scandata_field_list = ['/entry1/measurement', '/entry1/plotted']
#scan_command_field_list = ['/entry1/scan_command']
_entry = '/entry1'
_measurement = '/measurement'







class pdnx(pd.DataFrame): 
    '''
    nexusformat wrapper: tries to create dataframe from default data
    whole nexus file is under .nx attribute
    e.g. 
    n=pdnx(p % 633777)  open file for scan 633777 (p is filename/format specifier)
    n                   display pandas dataframe
    n.plot()            plot pandas dataframe
    n.plot('idgap','ic1monitor')	pandas plot with selected x and y collumns
    n.plt('idgap','ic1monitor')		same but with pdnx defaults (title etc)	
    n.nx                nexus tree
    print(n.nx.tree)     print nexus tree
    n.find('chi')	find 'chi' key(s) in tree and display value(s) (n.find() for all)
    n.findkeys('chi')	return list of key value lists for key 'chi'
    n.pruned_tree(n)    return nexus tree up to n levels deep
    n.nx.plot()         default nexus plot
    for i in range(633777, 633779):print(pdnx(p % i).scan)     print scan string for range of scans

    n['newkey'] = n.nx.entry1.before_scan.myval 	as long as 'newkey' is new then this pads out a new scan column with myval

    '''


    def __init__(self,  filestr, entry = _entry, data = _measurement):
        '''
        entry = select nexus entry for measurement data and set default to this entry
        data = nexus field containing datafor pandas dataframe
        '''
        try:
            _nx = nx.nxload(filestr,'r')

        except:
            print("=== Error loading file %s" % filestr)
            return
        
        _load_dataframe_success = False
        _use_classicscan = False

        entrydata = None
        if not entry == None and not data == None:
            entrydata = entry+data
        else:
            entry =  getNexusSubentryWithDefinition(_nx, definition = 'NXclassic_scan')
            for key in _nx[entry].keys():
                if 'NXdata' in str(type(_nx[entry][key])):
                    entrydata = '%s/%s' % (entry, key)
                    _use_classicscan = True
        print(entrydata)####################
        try:
            _nx[entrydata]
        except:
            raise ValueError('=== Problem finding data field')

        print(entrydata)


            
        try:
            if _use_classicscan:
                keys = _nx[entrydata]['scan_fields'] # use fields from scan_fields required by NXclassic_scan
            else:
                keys = _nx[entrydata].keys()        # use all fields - must all be the same length to avoid an error
            print(keys) #########################################
            #nx_scan_dict = dict(_nx[entry+data])
            nx_scan_dict = {}
            print(entrydata) #####################################
            #for key in nx_scan_dict.keys():
            for key in keys:
                print(key)################################
                try:
                    #nx_scan_dict[key] = nx_scan_dict[key].nxdata.flatten()
                    print(1, _nx[entrydata])
                    print(2, _nx[entrydata][key])
                    print(3, _nx[entrydata][key].nxdata.flatten())
                    print(4, len(_nx[entrydata][key].nxdata.flatten()))
                    
                    nx_scan_dict[key] = _nx[entrydata][key].nxdata.flatten() #### need ot uncomment
                    print('=== succeeded:', key)#########################
                except:
                    print('=== failed:', key)#########################
                    pass
            #print('dict---', nx_scan_dict)################### xx problem with dict - infinite regression
            pd.DataFrame.__init__(self, nx_scan_dict)   ###############
            _load_dataframe_success = True
        except:
            pass

        try:
            _nx['default'] = entry #set default entry to specified entry (for files with multiple enties)
        except:
            pass


        if not _load_dataframe_success:
            print('=== Failed to create DataFrame from data - create empty DataFrame')
            pd.DataFrame.__init__(self)

        setattr(self,'nx',_nx)
        

        try:
            setattr(self, 'scan', filestr+'\n'+_nx[entry]['title'].nxdata)
        except:
            pass

def getNexusSubentryWithDefinition(nxroot, definition = None):
    '''
    return NeXus tree branch string that is an entry or subentry containing the specified definition (string)
    if no definition specified then the function displays all the definitions found
    '''
    field_with_definition = None
    all_definitions = []

    for entry in nxroot.keys():                             # loop through each entry
        try:
            defn = str(nxroot[entry]['definition'])
            all_definitions += [defn]
            if defn == definition:
                #field_with_definition = nxroot[entry]           # if it has required definition then break
                field_with_definition = '/%s' % entry 
                break
        except:
            pass
    
        for subentry in nxroot[entry].keys():               # loop through each field in entry
            if 'NXsubentry' in str(type(nxroot[entry][subentry])): # if it is a subentry ...
                try:
                    defn = str(nxroot[entry][subentry]['definition'])
                    all_definitions += [defn]
                    if defn == definition:  # check if it has the required definition
                        #field_with_definition = nxroot[entry][subentry]
                        field_with_definition = '/%s/%s' % (entry, subentry) 
                        break
                except:
                    pass
                
    if definition == None:
        print(all_definitions)
        
    return field_with_definition


    def __init_old__(self,  filestr, entry = _entry, measurement = _measurement):
        '''
        entry = select nexus entry for measurement data and set default to this entry
        measurement = nexus field containing datafor pandas dataframe
        '''
        try:
            _nx = nx.nxload(filestr,'r')

        except:
            print("=== Error loading file %s" % filestr)
            return
        
        _load_dataframe_success = False

        try:
            nx_scan_dict = dict(_nx[entry+measurement])
            for key in nx_scan_dict.keys():
                try:
                    nx_scan_dict[key] = nx_scan_dict[key].nxdata.flatten()
                except:
                    pass
            pd.DataFrame.__init__(self, nx_scan_dict)
            _load_dataframe_success = True
        except:
            pass

        try:
            _nx['default'] = entry #set default entry to specified entry (for files with multiple enties)
        except:
            pass

        if not _load_dataframe_success:
            #print('=== Failed to create DataFrame from data - create empty DataFrame')
            pd.DataFrame.__init__(self)

        setattr(self,'nx',_nx)
        

        try:
            setattr(self, 'scan', filestr+'\n'+_nx[entry]['title'].nxdata)
        except:
            pass


    def plt(self, *args, **kwargs):
        kwargs.setdefault('title', self.scan)
        kwargs.setdefault('grid', True)
        self.plot(*args, **kwargs)


    def _list_to_dot_sep_string(self, lst):
        outstr = ''
        for item in lst:
            outstr += '.' + str(item)
        return outstr

    def _find_key(self, tree, key, previous_keys=[]):
        global _keylist
        try:
            for keyval in tree.keys():
                if keyval == key or key == '':
                    _keylist += [previous_keys + [keyval]]
                self._find_key(tree[keyval], key, previous_keys = previous_keys + [keyval])
        except:
            pass

    def findkeys(self, keystring):
        'Return list of key sequences (lists) that end with keystring'    
        global _keylist
        _keylist=[]
        self._find_key(self.nx, keystring)
        return _keylist

    def find(self, keystring=''):
        'Return nexus fields and values for keystring'
        for key_sequence in self.findkeys(keystring):
            obj = self.nx
            for key in key_sequence:
                obj = obj[key]
            print('.nx' + self._list_to_dot_sep_string(key_sequence) + ' : \t', obj)

    def pruned_tree(self, depth):
        'Print pruned tree'
        allfieldlist = self.findkeys('')
        previous = []
        for fieldlist in allfieldlist:
            fieldshort = fieldlist[:depth]
            if fieldshort != previous:
                print(self._list_to_dot_sep_string(fieldshort))
                previous = fieldshort





def vec2mat(vecx, vecy, vecz, n_inner=None):
    #matx, maty, matz = vec2mat(vecx, vecy, vecz, n_inner=None)
    #convert vectors from 2D scan to matrices
    #vecx,y,z: arrays (any dimension) or lists
    #matx,y,z: 2D arrays
    #n_inner: Number of points in inner loop - calculated if not specified
    #Arrays are truncated if the size doesn't match the required shape

    vx = np.array(vecx[:]); vy = np.array(vecy[:]); vz = np.array(vecz[:]) #get inputs in standard form
    if n_inner == None:   #calculate number in inner loop by looking for jumps
        jumps = np.abs(np.diff(vx) * np.diff(vy))
        n_inner = matplotlib.mlab.find(jumps>np.mean(jumps))[0] + 1

    n_outer = len(vx) // n_inner
    #reshape matrices
    matx = vx[0:n_inner * n_outer].reshape(n_outer,n_inner)
    maty = vy[0:n_inner * n_outer].reshape(n_outer,n_inner)
    matz = vz[0:n_inner * n_outer].reshape(n_outer,n_inner)

    return matx, maty, matz




############ testing - delete ##################
filenum = 815893
p='/dls/science/users/spc93/misc_nexus_data/modified/%i.nxs'
#n=pdnx(p % 815893)
#n=pdnx(p % 815893, entry = '/entry1/scan/', data = None)
n=pdnx(p % 815893, entry = None, data = None) # use classic_scan
#n=pdnx(p % 815893) # old method

##################

