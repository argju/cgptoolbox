"""
Basic example of cGP study using sigmoid model of gene regulatory networks.

This example is taken from Gjuvsland et al. (2011), the gene regulatory model 
is eq.(15) in that paper, and the example code reproduces figure 4a.

.. plot::
    :width: 400
    :include-source:

    from cgp.examples.sigmoid import cgpstudy, plot_gpmap
    gt, par, ph = cgpstudy()
    plot_gpmap(gt, ph)

The gene regulatory model is a
:class:`~cgp.sigmoidmodels.sigmoid_cvode.Sigmoidmodel`. It is a diploid model,
in the sense that there is two differential equations per gene, one for each
allele. This calls for a genotype-to-parameter map where the two alleles are
mapped to two independent sets of parameter values. The function
:func:`~cgp.gt.gt2par.geno2par_diploid` offers such genotype-to-parameter
mapping with allele-specific parameters .

:class:`~cgp.sigmoidmodels.sigmoid_cvode.Sigmoidmodel` encodes a system with
three genes and the parameter *adjmatrix* is used to specify connetivity and
mode of regulation. The model in Gjuvsland et al. (2011) only has two genes, so
we focus only on gene 1 and 2, and make sure that gene 3 is not regulating any
of them.

For a given parameter set the resulting phenotype is the equilibrium expression
level of gene 2. In order to find the equilibrium we combine
:class:`~cgp.sigmoidmodels.sigmoid_cvode.Sigmoidmodel` with
:class:`~cgp.phenotyping.attractor.AttractorMixin`.

Reference:

    * Gjuvsland AB, Vik JO, Wooliams JA, Omholt SW (2011)
      `Order-preserving principles underlying genotype-phenotype maps ensure
      high additive proportions of genetic variance <http://onlinelibrary.wiley.
      com/doi/10.1111/j.1420-9101.2011.02358.x/full>`_ Journal of Evolutionary
      Biology 24(10):2269-2279.
"""
# pylint: disable=W0621, W0212, W0612

import os
import numpy as np
import matplotlib.pyplot as plt
import cgp.sigmoidmodels

from cgp.utils.placevalue import Placevalue
from cgp.gt.gt2par import geno2par_haploid
from cgp.sigmoidmodels.sigmoid_cvode import Sigmoidmodel
from cgp.phenotyping.attractor import AttractorMixin
from cgp.utils.rec2dict import dict2rec

from collections import OrderedDict

########################################
### Model of gene regulatory network ###
########################################

class Model(Sigmoidmodel, AttractorMixin):
    """
    Class that combines a sigmoid model with the ability to find equilibria.
    
    .. inheritance-diagram:: Model cgp.sigmoidmodels.sigmoid_cvode.Sigmoidmodel cgp.phenotyping.attractor.AttractorMixin
       :parts: 1
    """
    pass

model = Model(adjmatrix=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]]))
model.par_mean[[0,6,12],:]=200 #single allele, increase lowest production

########################################
### Genotypes for two biallelic loci ###
########################################

names = ['Gene1', 'Gene2','Gene3']
genotypes = Placevalue([2, 2, 2], names=names, msd_first=False)

### Genotype-to-parameter map
gt2par = geno2par_haploid


#binary matrix indicating which parameter values are 
#inherited together with what locus
loc2par = np.zeros((3, 18)) #two loci, 18 parameters in total
loc2par[0, 0:6] = 1         
loc2par[1, 6:12] = 1
loc2par[2, 12:18] = 1

### sample parameter values and duplicate production levels
def sample_duplicatepar():
    par = model.sample_hetpar()
    p = par.view(dtype=np.float).reshape((18,2))
    p[:,1] = p[:,0]
    par.alpha1 = par.alpha1*[1,2]
    par.alpha2 = par.alpha2*[1,2]
    par.alpha3 = par.alpha3*[1,2]
    
    return par
    

### Parameter-to-phenotype map
def par2ph(par):
    """
    Phenotype: Equilibrium expression level of gene 2.
    
    This maps a parameter value to a phenotype value.
    """
    model._ReInit_if_required(y=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    model.pr[:] = par
    _t, y = model.eq(tmax=1000, tol=1e-6)
    return dict2rec(OrderedDict(zip(('Y11','Y12','Y21','Y22','Y31','Y32'),(y[0],y[1],y[2],y[3],y[4],y[5]))))

def get_ffl_phenotypes(par):
    motifnames = ['Lin1','Lin2','Lin3','Lin4','Coh1','Coh2','Incoh1','Incoh2']
    adjmatrix = [[0,0,0,1,0,0,0,1,0],
                 [0,0,0,1,0,0,0,-1,0],
                 [0,0,0,-1,0,0,0,-1,0],
                 [0,0,0,-1,0,0,0,1,0],
                 [0,0,0,1,0,0,1,1,0],
                      [0,0,0,-1,0,0,-1,1,0],
                      [0,0,0,1,0,0,1,-1,0],
                      [0,0,0,-1,0,0,1,-1,0]]
    #adjmatrix = np.genfromtxt('/home/gjuvslan/CompBio/cGPsandbox/sigmoidmodels/motifs_adjacency.csv',delimiter=',',skip_header=1)
    Y1 = np.nan*np.empty(len(adjmatrix))
    Y2 = np.nan*np.empty(len(adjmatrix))
    Y3 = np.nan*np.empty(len(adjmatrix))
    for i in range(len(Y3)):
        print motifnames[i],
        model.adjmatrix = np.array(adjmatrix[i]).reshape(3,3)
        model.boolfunc = model.adjacency_to_boolfunc(model.adjmatrix)
        P = par2ph(par)
        Y1[i] = P['Y11'][0]
        Y2[i] = P['Y21'][0]
        Y3[i] = P['Y31'][0]
    
    dtype = [(name, float) for name in motifnames]
     
    #return Y1.view(dtype, np.recarray), Y2.view(dtype, np.recarray), Y3.view(dtype, np.recarray) 
    return Y3.view(dtype, np.recarray)

def get_motif_adjacency_phenotypes(par):
    d = os.environ
    motifs = np.genfromtxt(os.path.dirname(cgp.sigmoidmodels.__file__)+'/motifs_adjacency.csv',delimiter=',',skip_header=1)
    subset = np.hstack((np.arange(motifs.shape[0],step=100),motifs.shape[0])) #divide motifs into subset of 100
    try:
        i = int(d['Nloci'])    #read in subset from bash variable (using the Nloci variable as a lazy hack)
    except KeyError:
        i = 1
    motifs = motifs[subset[i-1]:subset[i],]
    motifnames = range(subset[i-1]+1,subset[i]+1)    #motif number corresponding to line number in motifs_adjacency.csv (skip header)

    Y1 = np.nan*np.empty(motifs.shape[0])
    Y2 = np.nan*np.empty(motifs.shape[0])
    Y3 = np.nan*np.empty(motifs.shape[0])
    for i in range(len(Y3)):
        print i
        model.adjmatrix = np.array(motifs[i]).reshape(3,3)
        model.boolfunc = model.adjacency_to_boolfunc(model.adjmatrix)
        P = par2ph(par)
        Y1[i] = P['Y11'][0]
        Y2[i] = P['Y21'][0]
        Y3[i] = P['Y31'][0]
    
    dtype = [(str(name), float) for name in motifnames]
     
    #return Y1.view(dtype, np.recarray), Y2.view(dtype, np.recarray), Y3.view(dtype, np.recarray) 
    return Y3.view(dtype, np.recarray)


def cgpstudy():
    """
    Connecting the building blocks of a cGP study.
    
    This implements the simulation pipeline in Gjuvsland et al. (2011).
    """
    from numpy import concatenate as cat
    gt = np.array(genotypes)
    hetpar = sample_duplicatepar()
    par = cat([gt2par(list(g), hetpar, loc2par) for g in gt])
    #ph = [cat(y) for y in zip(*tuple([get_motif_adjacency_phenotypes(p) for p in par]))]
    #ph = cat([get_motif_adjacency_phenotypes(p) for p in par])
    ph = cat([get_ffl_phenotypes(p) for p in par])
    
    return gt, par, ph

if __name__ == "__main__":
    gt, par, ph = cgpstudy()
