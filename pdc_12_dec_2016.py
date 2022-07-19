from dlstools import dataloader
close('all')

path='/dls/i16/data/2016/cm14471-5/'
d=dataloader.dlsloader(path+'%i.dat')

d(619735)
print 'Time: %i:%i:%i' % (d.Hours[0], d.Minutes[0], d.Seconds[0])
'''
#
12/12/2016 PDC run

15:52 start data collection run 1: atten 255 eta 27.082 (etacen+delth)
17:38 end run 1 (beam lost) 

20:47 start run 2 (after waiting for beam) atten 220 eta 27.082 (etacen+delth)
23:53 end run 2

23:54 start run 3 atten 255 eta 27.107 (etacen+delth*1.5)
03:57 end run 3

03:58 start run 4 atten 220 eta 27.107 (etacen+delth*1.5)
08:03 end run 4

08:04 start run 5 (same conditions as run1 - data can be added to run 1)
08:33 end of run


#start 619734
#restart 619735 no beam for several hours between first two scans. script ended ok 8 a.m.
#restart script for last hour of beam 619739


pos atten 255; pos w 2; pos atten 255
pos delta 56
pos eta etacen+delth
scan x 1 1000 1 t 14 #3.8 hours

pos atten 220; pos w 2; pos atten 220
pos delta 56
pos eta etacen+delth
scan x 1 1000 1 t 14 #3.8 hours

pos atten 255; pos w 2; pos atten 255
pos delta 56
pos eta etacen+delth*1.5
scan x 1 1000 1 t 14 #3.8 hours

pos atten 220; pos w 2; pos atten 220
pos delta 56
pos eta etacen+delth*1.5
scan x 1 1000 1 t 14 #3.8 hours





'''

#savefig('/home/i16user/tmp/tmp.pdf')
