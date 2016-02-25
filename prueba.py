#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

import os

path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_PMW.fits"


fits_data = fits.getdata(path+filename,"image")
 
gauss = Gaussian2DKernel(stddev=2)
z = convolve(fits_data, gauss, boundary='extend')

pixels = None
fig, main_axes = plt.subplots()
main_axes.imshow(fits_data, origin="lower", interpolation="None")
#main_axes.imshow(z, origin="lower", interpolation="None")
plt.show()

#header = fits_data
#print header
	


