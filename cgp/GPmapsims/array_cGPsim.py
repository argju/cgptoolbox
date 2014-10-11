"""Python module designed for running cGP simulation studies as queue jobs."""
from GPutils import find_jobID, end_jobID
import numpy as np
import os
import pickle
import tables
import logging

# Logger objects keep a dict of some relevant information
keys = "asctime levelname name lineno process message"
logging.basicConfig(level=logging.INFO, format = "%(" + ")s\t%(".join(keys.split()) + ")s")
arraylog = logging.getLogger('array')

#import simulation settings from bash
d = os.environ
arraylog.info("** Task number: %s starting work" % d['TASK_ID'])

############################################
# Set up simulations for specific cGPmodel #
############################################

if d['Model']=='adjmotifshaploid':
    from cgp.examples.sigmoid_haploid import cgpstudy
    
else:
    Exception('Unknown cGPmodel: '+d['Model'])

########################################################
#  Case independent code for running array simulations #
########################################################

subdir = d.get('PBS_O_WORKDIR', d.get("SUBMITDIR")) #qsub or slurm queue system
os.chdir(subdir) #sbatch job starts in .tmp directory

jobID=0
while jobID<int(d['Nsims']):
    ## find next job
    jobID = find_jobID(d)
    
    ## do computations 
    if jobID<int(d['Nsims']):
        arraylog.info("** Starting solving rep nr: " + str(jobID))
        try:
            np.random.seed([int(d['Nloci']),jobID])
            #loci_names, hetpar, relvar = define_parameter_GPmap(int(d['Nloci']))
            #genotypes, parameters, phenotypes = cGPstudy(hetpar, relvar, geno2par, par2pheno, loci_names)
            genotypes, parameters, phenotypes = cgpstudy()
            with tables.openFile("%s_%s.hdf5" % (d["datafile"], jobID),'w') as f: #create hdf5 file for phenotype data 
                #f.createTable('/','hetpar', hetpar)
                #f.createArray('/','relvar', relvar)
                f.createTable('/', 'genotypes', genotypes)
                f.createTable('/', 'parameters', parameters)
                f.createTable('/', 'phenotypes', phenotypes)
            end_jobID(d,jobID)
        except:
            arraylog.exception('ERROR: exception caught during solution of replicate '+str(jobID))
            with open("%s_%s.pickle" % (d["errorfile"], jobID), 'w') as f:
                Nloci = d['Nloci']
                for i in Nloci, jobID:
                    pickle.dump(i, f)
        # remove timeout file
        arraylog.info("** Finished solving rep nr: " + str(jobID))
    
