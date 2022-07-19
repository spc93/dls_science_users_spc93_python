#	%run /dls_sw/i16/software/python/userscripts/steve/diffcalc_demo.py
#print 'hello'
# do this interactively...
module load diffcalc
diffcalc fourcircle


newub 'NaOsO3_Pnma'
setlat 'NaOsO3' 5.3927 7.6038 5.3507 90 90 90
ub
con a_eq_b
pos en 10.872
pos delta 30
pos eta 15
addref 1 1 1
showref
trialub  ???
setnhkl [0 0 1]
pos hkl [1 1 1]
checkub

setlat
#Bi
setlat 'Bi' 4.546 4.546 11.862 90 90 120

#006
pos delta 33.526
pos eta 1.6985
pos chi -35.435
pos phi -78.534
addref [0 0 6]

#02m2
pos delta 36.0
pos eta 8.0
pos chi -40.5
pos phi 96.2
addref [0 2 -2]

#012
pos delta 19.97
pos eta 8.65
pos chi -90
pos phi 0
addref [0 1 2]

