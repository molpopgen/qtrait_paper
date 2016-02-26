#This file defines some common functions for the following cases
# * simulations of "single genomic regions"

import fwdpy as fp,math

def make_neutral_region():
    return [fp.Region(0,1,1)]

def make_buried_rec_region(littler):
    """
    For a region of size littler in 
    a larger 50cM region, set up the 
    boundaries, assuming littler 
    corresponds to recombination
    along the interval [0,1)
    """
    rest = 0.5-littler
    ratio = rest/littler
    return {'region':[fp.Region(-ratio,1+ratio,1)],
            'beg':-ratio,
            'end':(1+ratio)}

def make_Gaussian_sregion(beg=0,end=1,weight=1,sigma = None):
    if sigma is None:
        raise RuntimeError("sigma cannot be None")
    return [fp.GaussianS(beg,end,weight,sigma)]

def get_sigE_additive(m,VS,H):
    ##Figure out sigma_E from params
    EVG = 4.0*m*VS*VS
    EVE = EVG*(1.0-H)/H
    sigE = math.sqrt(EVE)
    return sigE
