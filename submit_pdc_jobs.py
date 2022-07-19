import drmaa
import pickle
import time
from os import chmod
pickle_file='drmma_jobs_june_b_15.p'
#qsubopts='-h -l m_mem_free=20G -q low.q'  ########## -h= hold!
qsubopts='-l m_mem_free=20G -q low.q' #run on machine with 20Gb and ensure slow queue (no time limit) 
#qsubopts='-l m_mem_free=20G' # simplified options (
#to check jobs: 
#jobdict = pickle.load( open(pickle_file , "rb" ) )
#chkjobs(jobdict)
#
#Always use s for session; keep it open till a new one is needed

try: #close drmma session if open
    s.exit()
except:
    pass

start_time=time.time()

s = drmaa.Session(); s.initialize()
print 'A session was started successfully'
sessionstr = s.contact
print 'session contact returns: ' + sessionstr


def mk_txt_file(fullname, text):
    f = open(fullname, 'w')
    f.write(text)
    f.close()
    chmod(fullname,0755)

def chkjobs(sd):
    "dict: {'session': sessionstr, 'comlist': command string list, 'joblist': jobids}"
    #try: #close drmma session if open
    #    s.exit()
    #except:
    #    pass
    #s = drmaa.Session(sd['session']); s.initialize() #keep it open
    for ii in range(len(sd['joblist'])):
        print sd['session'], 'command= ', sd['comlist'][ii],'job_id= ', sd['joblist'][ii],'status=', s.jobStatus(jid)
    #s.exit()#keep it open
    print '%i jobs' % len(sd['joblist'])
    print 'Days since submitting jobs: %.2f' % ((time.time()-sd['start_time'])/3600/24)

pycommand='python'
pyfile='two_photon_correlations.py' #with extension
#pyargfmt='/dls/i16/data/2014/cm4968-3/ %i %i'
pyargfmt='/dls/i16/data/2014/cm4968-5/ %i %i'
bashfmt='submit_scan_%i.sh'
#path='/dls_sw/i16/software/gda/config/pythonscripts/data_analysis/src/'
path='/dls_sw/i16/software/python/userscripts/steve/'
datadir='/dls/i16/data/2014/cm4968-5/processing/'
bashdir='/home/spc93/tmp/bash/'
images_per_pt=9999 #2 for testing or 9999

#15/6/15
goodscns=[490213, 490216, 490219, 490222, 490226, 490230, 490231, 490232, 490790, 490793, 490796, 490799, 490802, 490805, 490808, 490811, 490814, 490817, 490820, 490823, 490826, 490829, 490832, 490835, 490838, 490841, 490844, 490847, 490850, 490853, 490856, 490859, 490862, 490865, 490868, 490871, 490874, 490877, 490880, 490883, 490886, 490889, 490892, 490895, 490898, 490901, 490904, 490907, 490910, 490913, 490916, 490919, 490922, 490925, 490928, 490931, 490934, 490937, 490940, 490943, 490946]
#goodscns=goodscns[0:1] #### quick test

jobids=[]
names=[]
for scan in goodscns:
    pyargs=pyargfmt % (scan, images_per_pt)
    bashname=bashfmt % scan
    print bashname
    fullname=bashdir+bashname
    #mk_txt_file(fullname,'module load scisoftpy/ana\ncd '+datadir+'\n'+pycommand+' '+path+pyfile+' '+pyargs)
    mk_txt_file(fullname,'. /etc/profile.d/modules.sh\n'+'module load scisoftpy/ana\ncd '+datadir+'\n'+pycommand+' '+path+pyfile+' '+pyargs) #test: add modules.sh
    print fullname
    print 'Creating job template'
    jt = s.createJobTemplate()
    jt.nativeSpecification=qsubopts
    jt.remoteCommand = fullname
    print 'run job'
    jid=s.runJob(jt)
    jobids+= [jid]
    names+=[fullname]
    s.deleteJobTemplate(jt)
#s.exit() #keep it open

jobdict={'session': sessionstr, 'comlist': names, 'joblist': jobids, 'start_time':start_time}
pickle.dump(jobdict, open(pickle_file, "wb" ) )
chkjobs(jobdict) # show job status (repeat for update)
