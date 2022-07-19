import sys
sys.path.append('/dls_sw/apps/diffcalc/3/diffcalc3')


from startup.i16 import * 

newub("Example_i16")
setlat("SiO2", 4.913, 5.405)
ub()
sixc()
pos(mu, -2)
addorient([0, 0, 1], [0, 0, 1], sixc(), "norm")
addorient([0, 1, 0], [0, 1, 0], sixc(), "plane")
setnhkl([1, 0, 0])
con(qaz, 0, alpha, 0, eta,  0)
hkl_to_angles(0, 0, 1)
angles_to_hkl((0.0, 0.0, 0.0, 3.31, 0.0, 10.62))
sim(hkl, [0, 0, 1])
pos(hkl, [0, 0, 1])
scan(alpha, 0, 10, 1, hkl, [0, 0, 1], sixc, hklverbose.alpha, hklverbose.beta)
scan(qaz, 0, 90, 10, hkl, [0, 0, 1], sixc, hklverbose.alpha, hklverbose.beta)
con(qaz, 90, psi, 90, chi, 90)
scan(psi, 90, 0, 10, hkl, [0, 0, 1], sixc, hklverbose.alpha, hklverbose.beta)
setmiscut(0)
scan(psi, 90, 0, 10, hkl, [0, 0, 1], sixc, hklverbose.alpha, hklverbose.beta)