"""
Module for handling job id's for embarrassingly parallel jobs 
(see array_cGPsim.py for use of this moduel). It handles reassignment of
timed out and error jobs.

Example:

Settings for the jobfile system are kept in a dict.
>>> d=dict(Nsims=3, jobfile='./job', lockfile='./lock', SimID='.', errorfile='./error')

The main functions are:
find_jobID : return next job ID to be done. If >= Nsims, all are done.
end_jobID : indicate that a job ID has finished.

find_jobID() returns the first unstarted job ID.
>>> jobID = find_jobID(d)
>>> jobID
0

end_jobID() removes the job's timeout file, indicating that it has finished.
>>> end_jobID(d,jobID)                            # job 0 completed

Get a new job ID.
>>> jobID = find_jobID(d)                      
>>> jobID
1

Failing to end the job indicates that it has timed out.
>>> # end_jobID(d,jobID)                          # job 1 timed out

Simulate an error by creating an "error" file.
>>> jobID = find_jobID(d)
>>> touch("error_2.pickle")                     # job 2 resulted in error file created in except statement
>>> end_jobID(d,jobID)                          # job 2 completed

# all jobs have been started once, start looking for timeouts and errors

Errored jobs...
>>> jobID = find_jobID(d)
>>> jobID
2
>>> end_jobID(d,jobID)                           # job 2 completed with success

...and timed-out jobs.
>>> jobID = find_jobID(d)
>>> jobID
1
>>> end_jobID(d,jobID)                           # job 1 completed with success

If find_jobID() returns a value >= Nsims, all jobs are done.
>>> jobID = find_jobID(d)
>>> jobID
3

Clean up after doctests. TODO: Use a temporary directory for testing.
>>> os.remove("job")
>>> os.remove("job_redone_2.pickle")

"""

from cgp.utils.poormanslock import Lock
import logging
import os
import shutil
import numpy as np

# Logger objects keep a dict of some relevant information
keys = "asctime levelname name lineno process message"
logging.basicConfig(level=logging.INFO, format = "%(" + ")s\t%(".join(keys.split()) + ")s")
arraylog = logging.getLogger('findjob')

def touch(path, times=None):
    """
    Portable alternative to Unix touch.
    
    >>> touch("__touchtest__")
    >>> os.remove("__touchtest__") # now the file exists
    """
    if not os.path.exists(path):
        open(path, "w").close()
    os.utime(path, times)

def find_jobID(d):
    """
    Find the next undone job. This function returns the next unfinished
    jobID such that 0 <= jobID < d['Nsims']. If all jobs have been assigned
    jobs ending with error or timed out jobs are redone. If all jobs are completed 
    without exception or timeout the the function returns d['Nsims']

    d['Nsims']     - total number of jobs
    d['jobfile']   - path to file containing the last started job_timeout
    d['lockfile']  - path to file containing the last started job_timeout
    d['SimID']     - simulation directory
    d['errorfile'] - path+filestem for errorfiles
    """
    
    #create jobfile if it does not exist
    if not os.path.exists(d['jobfile']):
        with open(d['jobfile'],"w") as f:
            f.write('0')
            
    #find a job id
    with Lock(lockname=d['lockfile'], retry_delay=0.2, max_wait=100):
        with open(d['jobfile']) as f:
            jobID = int(f.read())            

        if jobID<int(d['Nsims']):
            with open(d['jobfile'],"w") as f:
                touch("%s_timeout_%s" % (d["jobfile"], jobID))
                f.write(str(jobID+1))
            arraylog.info("Jobfile - next job: " + str(jobID))
            return jobID
        else:
            #redo jobs that exited with error or timed out
            files  = os.listdir(d['SimID'])
            np.random.shuffle(files)
            for file in files:
                if 'error' in file:
                    jobID = int(file.split('_')[-1].split('.')[0])
                    break # don't iterate over all the other files
            if jobID<int(d['Nsims']):
                touch("%s_timeout_%s" % (d["jobfile"], jobID))
                shutil.move("%s_%s.pickle" % (d["errorfile"], jobID), 
                            "%s_redone_%s.pickle" % (d["jobfile"], jobID))
                arraylog.info("Redoing failed job - next job: " + str(jobID))
                return jobID
            else:
                for file in files:
                    if 'job_timeout' in file:
                        jobID = int(file.split('_')[-1].split('.')[0])
            if jobID<int(d['Nsims']):
                touch("%s_timeout_%s" % (d["jobfile"], jobID))
                arraylog.info("Redoing timed out job - next job: %s" % jobID)
                return jobID

            else:
                arraylog.info("Jobfile - no jobs left - finishing ....")
                return d['Nsims']
                
def end_jobID(d,jobID):
    """
    Delete timeoutfile for jobID
    """
    os.remove("%s_timeout_%s" % (d["jobfile"], jobID))
       
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
