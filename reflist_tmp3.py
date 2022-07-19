import Crystal as Cr
cifpath='/home/spc93/spc_cifs/'

mc = Cr.Crystal()

#mc.load_cif(cifpath+'Bi2Se3_icsd_617072_cifbib.cif','617072-ICSD')
#mc.load_cif(cifpath+'Bi2Te3_icsd_658764_cifbib.cif','658764-ICSD')
#mc.load_cif(cifpath+'Al2O3_sapphire_icsd_10425_cifbib.cif','10425-ICSD')
#mc.load_cif(cifpath+'al2o3_sapphire_icsd_160607_cifbib.cif','160607-ICSD')
#mc.load_cif(cifpath+'Ba2CuHgO4.cif','1520216-COD')
mc.load_cif(cifpath+'LiNiPO4_icsd_184770_cifbib.cif','184770-ICSD')



refs=mc.reflection_list(8)
#eprev=0
#print; print '(h,k,l)\t\ttth\tIrel\td\tE(tth=180)'; print '-------\t\t---\t---'
#for ref in refs:
#    d=ref[5]
#    E180=12.4/2/d
#    if ref[2]>1.0 and eprev != E180:
#        print ref[0],'\t%.2f\t%.4g\t%.3f\t%.3f' % (ref[4], ref[2], d, E180)
#    eprev=E180
