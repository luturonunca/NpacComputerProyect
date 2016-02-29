#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import math
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel, convolve
import os
import library


plt.ion()

def main():
    ### SHORT WAVELENGH
    above_S, fits_data_S = library.get_data("S")
    image_S = above_S * fits_data_S
 
    # Create mask without NaN
    m_array_S = library.create_marray(fits_data_S)

            
    # Make the Mexican hat kernel and convolve
    # The values given to the kernels should be 1.7328, 2.3963, 3.5270
    # for S, M and L.
    mex_S = MexicanHat2DKernel(1.7328)
    mex_convol_S = convolve(image_S, mex_S, boundary = 'extend')
    m_mex_S = ma.masked_array(mex_convol_S, ~m_array_S.mask)

    ### MEDIUM WAVELENGH
    above_M, fits_data_M = library.get_data("M")
    image_M = above_M * fits_data_M
    
    # Create mask without NaN
    m_array_M = library.create_marray(fits_data_M)

            
    # Make the Mexican hat kernel and convolve
    # The values given to the kernels should be 1.7328, 2.3963, 3.5270
    # for S, M and L.
    mex_M = MexicanHat2DKernel(2.3963)
    mex_convol_M = convolve(image_M, mex_M, boundary = 'extend')
    m_mex_M = ma.masked_array(mex_convol_M, ~m_array_M.mask)





                        
    # Plot the figures; the min and the max values are reset in some of them
    # for visualization purposes.
    pixels = None
    fig, main_axes = plt.subplots(1,3, sharex=True, sharey=True, figsize=(15,5))
    #main_axes[0][0].imshow(mex, origin="lower", interpolation="None")
    #main_axes[0][0].set_title("MexHat kernel")
           
            
    main_axes[0].imshow(m_mex_S,origin="lower", interpolation="None")
    main_axes[0].set_title("Convolution with MexHat for S")
    

    main_axes[1].imshow(m_mex_M, origin="lower", interpolation="None")
    main_axes[1].set_title("Convolution with MexHat for M")
    
    #main_axes[1].imshow(m_gauss, origin="lower", interpolation="None")
    #main_axes[1].set_title("Convolution with gaussian")
           
           
    #main_axes[2].imshow(image, origin="lower", interpolation="None")
    #main_axes[2].set_title("Data")
            
    plt.show()
            

    #mean_data = (np.nanmax(fits_data) - np.nanmin(fits_data))/2
    #print "Mean", mean_data, "\nMax", np.nanmax(fits_data),\
    #      "\nMin", np.nanmin(fits_data)
    
    ## Plotting the histogram
    print "Histogram part"
    #mex_1d = m_mex.ravel()
    #histo = mex_1d
    #print histo2
    #histo = histo2[~np.isnan(histo2)]
    #indices = histo < 0
    #histo[indices] = 0
    #indices = histo > 0.004
    #histo[indices] = 0
    #print histo
    #bin_values, bin_boundaries = np.histogram(mex_1d[~np.isnan(mex_1d)], 200)
    #histo = m_mex.compressed()
    #print m_mex[1].compressed()
    #plt.hist(m_mex.compressed(), 200)
    #plt.show()


    #ma.compressed()
    #header = fits_data
    #print header
    return 0

if __name__ == '__main__':
    sys.exit(main())
