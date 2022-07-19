import TensorScatteringClass as ten
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.options.display.multi_sparse = False

minvalue = 1e-10


reflist = []
for h in range(6):
    for k in range(6):
        for l in range(6):
            reflist += [np.array([h, k, l])]
reflist.sort(key = lambda x : x[0]**2 + x[1]**2 + x[2]**2)


t=ten.TensorScatteringClass(CIFfile='ZnO Kisi et al icsd_67454.cif', Site='Zn1', verbose = False)


refdict = {}
for ref in reflist:

    t.TensorCalc(hkl = ref, K=2, Parity=1, Time=1)
    tensorAllowed = int(np.abs(np.dot(t.Fs, t.Fs)) > minvalue)
    t.gen_scalar_allowed
    t.site_scalar_allowed
    refdict[tuple(ref)] = {
            'SG_allowed': t.gen_scalar_allowed,
            'Site_allowed': t.site_scalar_allowed,
            'tensor_allowed': tensorAllowed
             }

df = pd.DataFrame(refdict).transpose()
df

