#!/usr/bin/env python
"""Reading in fits files"""

def loaddata(fitsfile):
    from astropy.io import fits
    hdulist = fits.open(fitsfile)
    
