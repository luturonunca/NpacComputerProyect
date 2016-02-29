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
from scipy.optimize import curve_fit

def gaussian_fit(w, max_value, mean_value, sigma):
    """
    Gaussian fit
    """
    # pylint: disable=E1101
    z = max_value*np.exp(-((w-mean_value)**2)/(2*(sigma**2)))
    return z
def get_gauss_sigma(data):
    bin_values, bin_boundaries = np.histogram(data, 200)
    my = np.float(np.max(np.nonzero(bin_values)))
    mx = np.float(np.max(bin_boundaries))
    normal_y = bin_values/my
    normal_x = bin_boundaries[:-1]/mx
    x = normal_x * mx
    y = normal_y * my
    fit, covariant = curve_fit(gaussian_fit, normal_x, normal_y)
    maxvalue = fit[0] * my
    background = fit[1] * mx
    dispersion = fit[2] * mx 
    return abs(dispersion)


def main():
    # Variables for the filename
    path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
    filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_PMW.fits"

    # Load the data
    fits_data = fits.getdata(path+filename,"coverage")
    # Remove zeros and NaNs from data
    data1d = np.ravel(fits_data)
    datanonan = data1d[~np.isnan(data1d)]
    data = datanonan[np.nonzero(datanonan)] 
    # Build histogram
    bin_values, bin_boundaries = np.histogram(data, 200)
    my = np.float(np.max(np.nonzero(bin_values)))
    normal_y = bin_values/my
    mx = np.float(np.max(bin_boundaries))
    normal_x = bin_boundaries[:-1]/mx
    x = normal_x * mx
    y = normal_y * my
 
    width = 0.7 * (bin_boundaries[1] - bin_boundaries[0])
    center = (bin_boundaries[:-1] + bin_boundaries[1:]) / 2
    
    fit, covariant = curve_fit(gaussian_fit, normal_x, normal_y)
    maxvalue = fit[0] * my
    background = fit[1] * mx
    dispersion = fit[2] * mx 
    sigma = dispersion 
    print "sigma: ", sigma," dispersion: ", dispersion 
    
    plt.plot(x, gaussian_fit(x, maxvalue, background, dispersion), 'r-', label='fit')
    
    plt.bar(center,bin_values , align='center', width=width) 
    yaux = [(-6. * dispersion) for i in range(len(x))]
    plt.plot(yaux,gaussian_fit(x, maxvalue, background, dispersion), 'g-') 

    print get_gauss_sigma(data)
    plt.show()


    return 0

if __name__ == '__main__':
    sys.exit(main())
