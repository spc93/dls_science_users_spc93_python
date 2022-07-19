
#import sys
#sys.path.insert(0,'/dls_sw/i16/software/python/crystal')

import Crystal as Cr

mycrys = Cr.Crystal()

mycrys.load_cif('/dls_sw/i16/software/python/crystal/CoO.cif','53058-ICSD')

mycrys.save_cif('/dls_sw/i16/software/python/crystal/test_save.cif')

test = Cr.Crystal()

test.load_cif('/dls_sw/i16/software/python/crystal/test_save.cif','53058-ICSD')

print test.general.cifstring

mycrys.reflection_list(5,refl='sym',anomalous_flag=False,sort='two_theta',print_list=True)

print '\n'

mycrys.save('/dls_sw/i16/software/python/crystal/test.json')

test2 = Cr.Crystal()

test2.load('/dls_sw/i16/software/python/crystal/test.json')

test2.reflection_list(5,refl='sym',anomalous_flag=False,sort='two_theta',print_list=True)





