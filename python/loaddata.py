#!/usr/bin/env python
"""Reading in fits files"""

def loaddata():
    from astropy.io import fits
    hdulist = fits.open('/Users/destry/Documents/Github/truffles_examples/GALFA_HI_RA+DEC_092.00+10.35_W.fits')
