import Crystal as Cr
cifpath='/home/spc93/spc_cifs/'

mcp = Cr.Crystal()
mcp.load_cif(cifpath+'FeBO3_icsd_34474.cif')# parent cif

mcm = Cr.Crystal()
mcm.load_cif(cifpath+'febo3_mgm3_p1.mcif')# magnetic cif



mcp.reflection_list(3, refl='allowed')
mcm.reflection_list(3, refl='allowed')







