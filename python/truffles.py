#!/usr/bin/env python  

"""3D machine Vision Source Finding"""

import numpy as np
from astropy.io import fits
from spore import spore

def truffles(FITS_FILENAME,verbose=False):
    if verbose == True:
        print 'VERBOSE!'
    hdulist = fits.open(FITS_FILENAME)
    cube = hdulist[0].data
    spore(cube, verbose=verbose)


#if __name__ == "__main__":
#    import sys
#    truffles(sys.argv[1])
