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

# Create masked array
def create_marray(image):
    mask = np.zeros((len(image), len(image[0])))
    for j in range (len(image)):
        for i in range (len(image[j])):
            if math.isnan(image[j][i]):
                continue
            else:
                mask[j][i] = 1

    # Create the masked array. 1-mask, be careful.
    m_array = ma.masked_array(image, 1-mask)
    return m_array

def main():
    
    # Variables for the filename
    path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
    filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_PMW.fits"
    
    # Load the data
    fits_data = fits.getdata(path+filename,"image")

    # Create masked array from fits_data
    m_array = create_marray(fits_data)
    
    # Make the Mexican hat kernel and convolve
    # The values given to the kernels should be 1.7328, 2.3963, 3.5270
    # for S, M and L.
    mex = MexicanHat2DKernel(2.3963)
    w = convolve(fits_data, mex, boundary = 'extend')
    m_mex = ma.masked_array(w, m_array.mask)
    #c_mex = np.multiply(w, m_array.mask)
    
    # Make the gaussian kernel and convolve
    gauss = Gaussian2DKernel(stddev=3)
    z = convolve(fits_data, gauss, boundary='extend')
    m_gauss= ma.masked_array(w, m_array.mask)
    #c_gauss = np.multiply(z, m_array.mask)
    
    # Plot the figures; the min and the max values are reset in some of them
    # for visualization purposes.
    #pixels = None
    #fig, main_axes = plt.subplots(2,2)
    #main_axes[0][0].imshow(mex, origin="lower", interpolation="None")
    #main_axes[0][0].set_title("MexHat kernel")
   
    
    #main_axes[0][1].imshow(m_mex, vmin=-0.1, vmax=0.05,\
    #                       origin="lower", interpolation="None")
    #main_axes[0][1].set_title("Convolution with MexHat")
    
    #main_axes[1][0].imshow(m_gauss, origin="lower", interpolation="None")
    #main_axes[1][0].set_title("Convolution with gaussian")
   
   
    #main_axes[1][1].imshow(fits_data, origin="lower", interpolation="None")
    #main_axes[1][1].set_title("Data")
    
    #plt.show()
    
    mean_data = (np.nanmax(fits_data) - np.nanmin(fits_data))/2
    print "Mean", mean_data, "\nMax", np.nanmax(fits_data),\
          "\nMin", np.nanmin(fits_data)
    
    ## Plotting the histogram
    print "Histogram part"
    mex_1d = m_mex.ravel()
    histo = mex_1d[~np.isnan(mex_1d)]
    indices = histo < 0
    histo[indices] = 0
    indices = histo > 0.004
    histo[indices] = 0
    print histo
    #bin_values, bin_boundaries = np.histogram(mex_1d[~np.isnan(mex_1d)], 200)
    plt.hist(histo, 200)
    plt.show()


    #ma.compressed()
    #header = fits_data
    #print header
    return 0

if __name__ == '__main__':
    sys.exit(main())
