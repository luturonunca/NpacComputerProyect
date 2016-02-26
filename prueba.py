#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel, convolve
import os

# Variables for the filename
path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_PMW.fits"

# Load the data
fits_data = fits.getdata(path+filename,"image")

# Make the Mexican hat filter
mex = MexicanHat2DKernel(2, x_size=5, y_size=5)
z = convolve(fits_data, mex, boundary = 'None')

# Make the gaussian filter
#gauss = Gaussian2DKernel(stddev=2)
#z = convolve(fits_data, gauss, boundary='extend')

# Plot the figure
pixels = None
fig, main_axes = plt.subplots()
#main_axes.imshow(fits_data, origin="lower", interpolation="None")
main_axes.imshow(z, origin="lower", interpolation="None")
plt.show()

#header = fits_data
#print header
