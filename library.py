#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
personal library
"""
import sys                                                                      
import math                                                                     
import numpy as np                                                              
import numpy.ma as ma                                                           
import matplotlib.pyplot as plt                                                 
from astropy.io import fits
from astropy.wcs import WCS                                                     
import os
from scipy.optimize import curve_fit


def transform(x, y):
  
    with fits.open('specific.fits') as fits_data:
        # get header
        header = fits_data[0].header
        wcs = library.WCS(header)
        ra_x , dec_y = wcs.convert_to_radec(x, y)
        ra, dec = '%.3f' % ra_x, '%.3f' % dec_y
    return ra, dec

def gaussian_fit(w, max_value, mean_value, sigma):
    """
    Gaussian fit
    """
    # pylint: disable=E1101
    z = max_value*np.exp(-((w-mean_value)**2)/(2*(sigma**2)))
    return z


def get_gauss_sigma(data):
    """
    returns the sigma of a gaussian fit over the fluxs
    """
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


def create_marray(image):
    """
    gets data to mask
    creates mask without NaN
    returns mask
    """
    mask = np.zeros((len(image), len(image[0])), dtype=bool)
    for j in range (len(image)):
        for i in range (len(image[j])):
            if math.isnan(image[j][i]):
                continue
            else:
                mask[j][i] = True

    # Create the masked array. 1-mask, be careful.
    m_array = ma.masked_array(image, mask)
    return m_array


def get_data(W):
    """
    gets either S, M or L as string representing the dataset wavelengh
    returns: the cutted data, i.e centered, and a boolean matrix
    that represent a threathole in the coverage matrix.
    
    """
    
    path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
    filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_P"\
    +W+"W.fits"
    # Load the data
    fits_data_uncut = fits.getdata(path+filename,"image")
    
    fits_coverage_uncut = fits.getdata(path+filename,"coverage")
    # Define cuts
    ycenter = int(len(fits_data_uncut)/2)
    xcenter = int(len(fits_data_uncut[0])/2)
    cut = int(ycenter/2)
    print "centers ", float(xcenter)/ycenter
    # Cut array to focus on center
    fits_data = fits_data_uncut[xcenter-cut:xcenter+cut,ycenter-cut:ycenter+cut]
    fits_coverage = fits_coverage_uncut[xcenter-cut:xcenter+cut,ycenter-cut:ycenter+cut]

    # Create masked array from fits_data
    m_array = create_marray(fits_data)
    # Create masked array from fits_data maskin with coverage
    mask_coverage = m_array.mask * fits_coverage
    above = fits_coverage > 8

    image = (above * fits_data)
    
    return above, fits_data



def transform(array,W):
    path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
    filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_P"\
    +W+"W.fits"
    header = fits.getheader(path+filename)
    w = WCS(header)
    array2 = w.wcs_pix2world(array,0)
    return array2


