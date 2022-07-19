import Crystal as Cr
cifpath='/home/spc93/spc_cifs/'

mc = Cr.Crystal()

#mc.load_cif(cifpath+'Bi2Se3_icsd_617072_cifbib.cif','617072-ICSD')
#mc.load_cif(cifpath+'Bi2Te3_icsd_658764_cifbib.cif','658764-ICSD')
#mc.load_cif(cifpath+'Al2O3_sapphire_icsd_10425_cifbib.cif','10425-ICSD')
mc.load_cif(cifpath+'MoS2_icsd_49801_cifbib.cif','49801-ICSD')

energy=12.4/1.54
refs=mc.reflection_list(energy)

print; print '(h,k,l)\t\ttth\tIrel'; print '-------\t\t---\t---'
for ref in refs:
    print ref[0],'\t%.2f\t%.2f' % (ref[4], ref[2])

#(1,0,l)
#look at 10l 20l 11l
shortlist=[ref for ref in refs if ref[0][0]==0 and ref[0][1]==0]
print; print '(h,k,l)\t\ttth\tIrel'; print '-------\t\t---\t---'
for ref in shortlist[0:25]:
    print ref[0],'\t%.2f\t%.2f' % (ref[4], ref[2])


#refs=mc.reflection_list(energy)
#refs=[ref for ref in refs if ref[0][0]==ref[0][1]==ref[0][2]] #filter refs hhh
print; print '(h,k,l)\t\ttth\tIrel\td\tE(tth=90)'; print '-------\t\t---\t---'
for ref in shortlist[0:25]:
    d=ref[5]
    E90=12.4/sqrt(2)/d
    print ref[0],'\t%.2f\t%.4g\t%.3f\t%.3f' % (ref[4], ref[2], d, E90)
