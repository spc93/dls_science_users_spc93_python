from tripod_test_new import *
from time import sleep

# Example script for tripod testing
# type  '%run tripod_test_demo' in ipython to run
# commands can be typed from command line (then print is not needed)

# show()			display all tripod settings
# mv_angs([ang1, ang2, ang3])	move tripod to these angles keeping translations fixed
# mv_trans([x,y,z])		move tripod to these tool point coordinates keeping angles fixed
# ctool()			return toolpoint xyz and angles
# q				stop all motors (no brackets but need to hit return)
# a				stop all motors (no brackets but need to hit return)
# z				stop all motors (no brackets but need to hit return)

# show settings
print show()
sleep(5)








'''
# move third angle to 0.5 deg, others zero
mv_angs([0,0,.5])
print show()
sleep(5)

# move z by 1 mm
C, alpha = ctool()	#get current values
mv_trans([C[0], C[1], C[2]+1.0])	# do the move
print show()		#show new settings
sleep(5)

# move z by -1 mm
C, alpha = ctool()	#get current values
mv_trans([C[0], C[1], C[2]-1.0])	# do the move
print show()		#show new settings
sleep(5)

# move angles back to zero
mv_angs([0,0,0])
print show()

'''


#Testing VFM pitch stability 19 June 2015 

# move angles to zero
mv_angs([0,0,0])
print show()
sleep(5)



# move x by 1 mm
C, alpha = ctool()	#get current values
mv_trans([C[0]+2, C[1], C[2]])	# do the move
print show()		#show new settings
sleep(5)

# move x by -1 mm
C, alpha = ctool()	#get current values
mv_trans([C[0]-2, C[1], C[2]])	# do the move
print show()		#show new settings
sleep(5)

