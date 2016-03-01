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
# Create masked array

def main():
    wavelength_range = sys.argv[1]
    print wavelength_range
    above, fits_data = library.get_data(wavelength_range)
    image = above * fits_data

    # Create mask without NaN
    m_array = library.create_marray(fits_data)
         
    # Make the Mexican hat kernel and convolve
    # The values given to the kernels should be 1.7328, 2.3963, 3.5270
    # for S, M and L.
    mex = MexicanHat2DKernel(1.7328)
    mex_convol = convolve(image, mex, boundary = 'extend')
    m_mex = ma.masked_array(mex_convol, ~m_array.mask)
    #c_mex = np.multiply(w, m_array.mask)
            
    # Make the gaussian kernel and convolve
    gauss = Gaussian2DKernel(stddev=1.7328)
    gauss_convol = convolve(image, gauss, boundary='extend')
    m_gauss= ma.masked_array(gauss_convol, ~m_array.mask)
    #c_gauss = np.multiply(z, m_array.mask)
            
    # Plot the figures; the min and the max values are reset in some of them
    # for visualization purposes.
    pixels = None
    fig, main_axes = plt.subplots(1,3, sharex=True, sharey=True, figsize=(15,5))
    #main_axes[0][0].imshow(mex, origin="lower", interpolation="None")
    #main_axes[0][0].set_title("MexHat kernel")
           
            
    main_axes[0].imshow(m_mex,\
                           origin="lower", interpolation="None")
    main_axes[0].set_title("Convolution with MexHat")
            
    main_axes[1].imshow(m_gauss, origin="lower", interpolation="None")
    main_axes[1].set_title("Convolution with gaussian")
           
           
    main_axes[2].imshow(image, origin="lower", interpolation="None")
    main_axes[2].set_title("Data")
            
    plt.show()
            

    mean_data = (np.nanmax(fits_data) - np.nanmin(fits_data))/2
    print "Mean", mean_data, "\nMax", np.nanmax(fits_data),\
          "\nMin", np.nanmin(fits_data)
    
    


    #ma.compressed()
    #header = fits_data
    #print header
    return 0

if __name__ == '__main__':
    sys.exit(main())
