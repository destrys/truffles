#!/usr/bin/env python                                                                                         
"""Reading in fits files"""

def example1():

    import sys
    import scipy as sp
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import pyfits

   #Path to fits file to be imported                                                                         
    data1 = "/Users/destry/Documents/Github/truffles_examples/GALFA_HI_RA+DEC_092.00+10.35_W.fits"

    #Read out basic info                                                                                      
    pyfits.info(data1)

    #Load file header into keys                                                                               
    header = pyfits.getheader(data1)
    header.keys()

    #Load actual data                                                                                         
    data_cube = pyfits.getdata(data1, 0)

    print 'Type: ', type(data_cube)
    print 'Shape:', data_cube.shape

    #I'm just going to look at a random slice                                                                 
    slice1 = data_cube[45, :, :]

    #Show the slice                                                                                           
    plt.imshow(slice1)
    plt.winter()
    plt.show()

