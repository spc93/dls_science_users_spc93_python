#use V-F on diode
#checkbeam and V-F?????
# measure a couple of Bragg peaks to check OK and set slit width 

phi_nosample=0
sx_nosample=0
sy_nosample=0
tsec=300

refs=[[1,1,1],[2,2,0],[3,1,1],[2,2,2],[4,0,0],[3,3,1],[4,2,2],[5,1,1],[3,3,3],[4,4,0],[5,3,1],[4,4,2],[6,2,0],[5,3,3],[6,2,2],[4,4,4],[7,1,1],[5,5,1],[6,4,2],[7,3,1],[5,5,3],[8,0,0],[7,3,3],[6,4,4],[8,2,2],[6,6,0]]

a_diamond=3.56699
a_silicon=5.43053

def calc_d(a,hkl):
    return a/sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2)

#test scans
#pos energy 12
#scancn idgap 0.002 21 w .5 ic1; go maxval
#pos delta 2*180/pi*arcsin(6.2/en()/calc_d(a_diamond,[4,4,0]))
#scancn delta .005 21 21 t .5


param_list=[{'name':'diamond_6um','a':a_diamond, 'sx':0.1, 'sy':-0.2, 'phi': 20.0, 'refs':refs},
        {'name':'diamond_3um','a':a_diamond, 'sx':0.1, 'sy':-0.2, 'phi': 20.0,'refs':refs},
        {'name':'silicon_nbs640','a':a_silicon, 'sx':0.3, 'sy':0.2, 'phi': 20.0, 'refs':refs}]

ntot=len(refs)*len(param_list)
print "Total number of reflections measured: %i" % ntot 
print "Seconds per reflection: %.0f" % tsec
print "Total hours: %.1f" % (tsec*ntot/3600.0)


crash

tthval=120.0 #xxxxx as large as possible
emax=15.0
emin=4.8
#pos uharmonic xxxxx


for params in param_list:   #loop through samples
    print params

    for ref in params['refs']:      #loop through reflections from each sample
        d=calc_d(params['a'],ref)
        e_tth=6.2/d/sin(tthval/2*pi/180)
        e_back=6.2/d
        if e_back>emin and e_tth<emax:  #measure if energy OK for Bragg peak and Bragg edge
            print ref, e_back, e_tth
            
            print 'Bragg peak scans'
            #pos sx params['sx'] sy params['sy'] phi params['phi']
            #pos energy e_tth
            #note('idgap scan')
            #scancn idgap 0.002 21 w .5 ic1; go maxval
            tag=params['name']+' '+str(ref)+' '+'Bragg peak'; print tag; note(tag);
            #pos delta tthval tthp tthp.apd
            #scancn energy 0.002 21 t .5
            
            print 'Bragg edge scans'
            #pos energy e_back sx params['sx'] sy params['sy'] phi params['phi']
            tag='idgap scan'
            #scancn idgap 0.002 21 w .5 ic1; go maxval
            tag=params['name']+' '+str(ref)+' '+'Bragg edge'; print tag; note(tag);
            #pos delta 0 tthp tthp.diode
            #scancn energy 0.002 21 t .5 #make longer and finer
            tag=params['name']+' '+str(ref)+' '+'no sample'; print tag; note(tag);
            #pos sx sx_nosample sy sy_nosample phi phi_nosample
            #scancn energy 0.002 21 t .5 #make longer and finer
