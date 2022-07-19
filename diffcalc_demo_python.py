#execfile('/dls_sw/i16/software/python/userscripts/steve/diffcalc_demo_python.py')
#diffcalc --python fourcircle

print 'hello'

newub('Bi')
while True:
	try:
		delref(1)
	except:
		break
		
setlat('Bi', 4.546, 4.546, 11.862, 90, 90, 120)

thoffset=2.2
#thoffset=-3

pos(en, 10.872)
#006
pos(delta, 33.526)
pos(eta, 1.6985+thoffset)
pos(chi, -35.435)
pos(phi, -78.534)
addref([0, 0, 6])

#02m2
pos(delta, 36.0)
pos(eta, 8.0+thoffset)
pos(chi, -40.5)
pos(phi, 96.2)
addref([0, 2, -2])

#012
pos(delta, 19.97)
pos(eta, 8.65+thoffset)
pos(chi, -90)
pos(phi, 0)
addref([0, 1, 2])

#swapref(2,3)

ub
#showref()
checkub()


