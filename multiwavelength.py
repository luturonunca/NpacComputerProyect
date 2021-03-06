#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import math
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel, convolve
from astropy.coordinates import Angle
from wcsaxes import WCS
import os
import library


plt.ion()


def get_filename(w):
    path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
    filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_P"\
    +w+"W.fits"
    return path+filename

def main():
    fig = plt.figure()




     

    ### LARGE WAVELENGH
    above_L, fits_data_L = library.get_data("L")
    image_L = above_L * fits_data_L
    _ , header_L = fits.getdata(get_filename("L"),header=True)

    # Ploting section
    wcs_proj_L = WCS(header_L)
    ax_wcs_L = fig.add_subplot(1,3,3, projection = wcs_proj_L)
    # stablish limits of plots
    # that work for the other plots
    world = wcs_proj_L.wcs_pix2world(57.,57., 0)
    zero = wcs_proj_L.wcs_pix2world(0.,0., 0)
    # build limits
    pix_L = wcs_proj_L.wcs_world2pix(world[0],world[1], 0)
    pix0_L = wcs_proj_L.wcs_world2pix(zero[0],zero[1], 0)
    # set limits
    ax_wcs_L.set_xlim(pix0_L[0], pix_L[0])
    ax_wcs_L.set_ylim(pix0_L[1], pix_L[1])
    # plots in sky coordinates
    ax_wcs_L.imshow(image_L, origin='lower', interpolation='None')
    ax_wcs_L.coords['ra'].set_ticks(color='red')
    ax_wcs_L.coords['dec'].set_ticks(color='red')


    ### SHORT WAVELENGH
    above_S, fits_data_S = library.get_data("S")
    image_S = above_S * fits_data_S
    _ , header_S = fits.getdata(get_filename("S"),header=True)
    # Create mask without NaN
    m_array_S = library.create_marray(fits_data_S)
   
    # check limits 289.345018791 -33.6099453587
    # Ploting section
    wcs_proj = WCS(header_S)
    ax_wcs = fig.add_subplot(1,3,1, projection = wcs_proj) #### OJO!!
    # defines coordinates of limits
    # -flag- print "world ", world[0], world[1]
    pix = wcs_proj.wcs_world2pix(world[0],world[1], 0)
    pix0 = wcs_proj.wcs_world2pix(zero[0],zero[1], 0)  
    #print "pix ", pix[0], pix[1] 
    ax_wcs.set_xlim(pix0[0], pix[0])
    ax_wcs.set_ylim(pix0[1], pix[1])
    ax_wcs.imshow(image_S, origin='lower', interpolation='None')
    ax_wcs.coords['ra'].set_ticks(color='red')
    ax_wcs.coords['dec'].set_ticks(color='red')

    ### MEDIUM WAVELENGH
    above_M, fits_data_M = library.get_data("M")
    image_M = above_M * fits_data_M
    _ , header_M = fits.getdata(get_filename("M"),header=True)
    # Create mask without NaN
    m_array_M = library.create_marray(fits_data_M)
     
    # Ploting section
    wcs_proj_M = WCS(header_M)
    ax_wcs_M = fig.add_subplot(1,3,2, projection = wcs_proj_M)
    # define limits
    pix_M = wcs_proj_M.wcs_world2pix(world[0],world[1], 0)
    pix0_M = wcs_proj_M.wcs_world2pix(zero[0],zero[1], 0)
    #print "pix ", pix[0], pix[1] 
    ax_wcs_M.set_xlim(pix0_M[0], pix_M[0])
    ax_wcs_M.set_ylim(pix0_M[1], pix_M[1])
    ax_wcs_M.imshow(image_M, origin='lower', interpolation='None')
    ax_wcs_M.coords['ra'].set_ticks(color='red')
    ax_wcs_M.coords['dec'].set_ticks(color='red')
    pix = wcs_proj_M.wcs_world2pix(world[0],world[1], 0)
    # -flag- print "pix ", pix[0], pix[1] 
    
    plt.show()

    return 0

if __name__ == '__main__':
    sys.exit(main())
