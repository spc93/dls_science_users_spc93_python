from dlstools import dataloader
close('all')

path='/dls/i16/data/2017/mt16226-1/'
d=dataloader.dlsloader(path+'%i.dat')

for scan in range(637412,637417+1): 
    d(scan)
    d.plot('ic1monitor');  

#savefig('/home/i16user/tmp/tmp.pdf')


'''
eta_offset=0.5
start:    18:23 5/5/17    Beam off for ~ 1 hour
end:    06:25    6/5/17
start:    06:31 7/5/17
end:    18:32 7/5/17


eta_offset=0.25
start:    06:26 6/5/17    Beam off for ~ 1 hour (at end)
end:    18:27    6/5/17
start: 18:33 7/5/17
end:    06:34 8/5/17

eta_offset=1.0
start:    18:28    6/5/17    Beam off for ~ 1 hour (at start)
end:    06:30    7/5/17
start: 06:35 8/5/17
end: 08:54    8/5/17 (Beam off at end)
    
    
'''