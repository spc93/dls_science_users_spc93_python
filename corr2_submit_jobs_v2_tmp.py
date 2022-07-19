# test to see when jobs go from queue and test to see if .e file is empty
# update dict with completion status: completed; success
# submit all jobs with -h option (change two_photon_... from test mode)

#drmaa - python job management
#qacct -j 3778077 #check job info



import os
from subprocess import Popen, PIPE
import pickle
import time
import pprint

MonitorOnly=False    #Monitor jobs already submitted - no new jobs

def oscall(cmd, test=False):
    print cmd
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    return p.communicate()

def mk_txt_file(fullname, text):
    f = open(fullname, 'w')
    f.write(text)
    f.close()

pycommand='python'
pyfile='two_photon_correlations_v2.py' #with extension
pyargfmt='/dls/i16/data/2014/cm4968-3/ %i'
#qsubopts='-pe smp 12' #for more cores (=more memory)
qsubopts='-l m_mem_free=20G -q low.q' #run on machine with 20Gb and ensure slow queue (no time limit) 
#qsubopts='-h -l m_mem_free=20G -q low.q' # hold till released
success_message='=== Job completed ==='
#pickle_file='jobs_aug_14.p' #for jobs dict
pickle_file='jobs_aug_14.p' #

bashfmt='submit_scan_%i.sh'

path='/dls_sw/i16/software/gda/config/pythonscripts/data_analysis/src/'
datadir='/dls/i16/data/2014/cm4968-5/processing/'

goodscns=[490213, 490216, 490219, 490222, 490226, 490230, 490231, 490232, 490790, 490793, 490796, 490799, 490802, 490805, 490808, 490811, 490814, 490817, 490820, 490823, 490826, 490829, 490832, 490835, 490838, 490841, 490844, 490847, 490850, 490853, 490856, 490859, 490862, 490865, 490868, 490871, 490874, 490877, 490880, 490883, 490886, 490889, 490892, 490895, 490898, 490901, 490904, 490907, 490910, 490913, 490916, 490919, 490922, 490925, 490928, 490931, 490934, 490937, 490940, 490943, 490946]
goodscns=goodscns[0:2]

if MonitorOnly==True:
    jobs = pickle.load( open(pickle_file , "rb" ) )
else:
    jobs={}
    #for scan in goodscns[0:1]:
    for scan in goodscns:
        pyargs=pyargfmt % scan
        bashname=bashfmt % scan
        mk_txt_file(datadir+bashname, 'module load scisoftpy/ana\ncd '+datadir+'\n'+pycommand+' '+path+pyfile+' '+pyargs)
        #os.system('cd '+datadir+'; qsub '+qsubopts+' '+datadir+bashname)
        cmd='cd '+datadir+'; qsub '+qsubopts+' '+datadir+bashname
        stdout, stderr = oscall(cmd)
        job_no=int([token for token in stdout.split() if token.isdigit()][0])
        jobs[job_no]= {'scan':scan, 'stdout':stdout, 'stderr':stderr}
    pickle.dump(jobs, open( "jobs_aug_14.p", "wb" ) )

success=0
for key in jobs.keys():
    scan = jobs[key]['scan']
    job_out_file = (datadir+bashfmt+'.o%i') % (scan, key)
    print job_out_file
    jobs[key]['success']=None
    try:
        f=open(job_out_file); lines=f.readlines(); f.close;
        for line in lines:
            if success_message in line:
                jobs[key]['success']=True
                success+=1
                pickle.dump(jobs, open( "jobs_aug_14.p", "wb" ) )
    except:
        print '=== Error reading file: '+job_out_file
    
pprint.pprint(jobs)
strlist=oscall('qstat')[0].split('\n')[2:-1]
livejobs=[int(str.split()[0]) for str in strlist]
print '=== Total jobs: %i \t Live jobs: %i \t Completed: %i' % (len(jobs), len(livejobs), success)

