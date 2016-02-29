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





