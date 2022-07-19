import Crystal as Cr
cifpath='/home/spc93/spc_cifs/'

mc = Cr.Crystal()
#mc.load_cif(cifpath+'CdTe_F43m_icsd_620529.cif','620529-ICSD')
#mc.load_cif(cifpath+'K_Br_Fm3m_icsd_44282.cif','44282-ICSD')
#mc.load_cif(cifpath+'Mo_Im3m_icsd_643957.cif','643957-ICSD') #004 good for Ir L3 11.2 keV  - have Mo 002
mc.load_cif(cifpath+'W_Im3m_icsd_653430.cif','653430-ICSD')#004 good for Ir L3 11.2 keV

ll=mc.reflection_list(15)

print; print 'Energy for 45 degree reflection (keV)'
for ref in ll:
    en=12.4/sqrt(2)/ref[5]
    print ref[0],'\t', en







