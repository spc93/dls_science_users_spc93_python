import drmaa

s = drmaa.Session(); s.initialize()
print 'A session was started successfully'
response = s.contact
print 'session contact returns: ' + response

#fullname='/dls/i16/data/2014/cm4968-3/processing/test_script3.sh'
fullname='/dls/i16/data/2014/cm4968-3/processing/submit_scan_458475.sh'

print 'Creating job template'
jt = s.createJobTemplate()
jt.remoteCommand = fullname

jobid = s.runJob(jt)
print 'Your job has been submitted with id ' + jobid

print '=== Job status', s.jobStatus(jobid)

#retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
#print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited)
#print retval


print 'Cleaning up'
s.deleteJobTemplate(jt)


s.exit()
print 'Exited from session'