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

def get_FWHM(wl):
    if wl == "S":
        fwhm = 17.6
    elif wl == "M":
        fwhm = 23.9
    elif wl == "L":
        fwhm = 35.2
    else:
        print "ERROR!! not a know image"
        fwhm = -1000
    return fwhm

def get_pixFWHM(wl):
    """
    gets Full Width Maximun for the diferent images
    """
    if wl == "S":
        fwhm = 1.7328
    elif wl == "M":
        fwhm = 2.3963
    elif wl == "L":
        fwhm = 3.5270
    else:
        print "ERROR!! not a know image"
        fwhm = -1000
    return fwhm

def peakvsback(data,centroid,wl,header):
    """
    makes de integral of the galaxy centered at 
    centroid
    """   
    wcs = WCS(header)
    fwhm = get_pixFWHM(wl)
    cent_pix = wcs.wcs_world2pix(centroid[0], centroid[1],0) 
    r_ap = 2.5 * fwhm 
    r_int = 3.5 * fwhm
    r_ext = 5.5 * fwhm
    # set iteration limits
    # x_l = left; x_r = right
    # y_u = up; y_b = bottom
    x_l = int(cent_pix[0] - r_ext)
    x_r = int(cent_pix[0] + r_ext)
    y_u = int(cent_pix[1] + r_ext)
    y_b = int(cent_pix[1] - r_ext)
    integral = 0
    background = 0
    npix = 0
    for j in range(y_b, y_u):
        for i in range(x_l, x_r):
            x_local = i - cent_pix[0]
            y_local = j - cent_pix[1] 
            r = np.sqrt((x_local**2)+(y_local**2))
            if r <= r_ap:
                integral +=  data[j][i]
            elif r >= r_int and r <= r_ext:
                npix += 1.
                background += data[j][i]  
    return integral, (background/npix)


def main():
    fig = plt.figure()


    ### LARGE WAVELENGH
    above_L, fits_data_L = library.get_data("L")
    image_L = above_L * fits_data_L
    _ , header_L = fits.getdata(get_filename("L"),header=True)
    cento = [289.3972219, -33.60861111]
    coco = peakvsback(image_L, cento, "L" , header_L)
    print " coco ", coco


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
    ax_wcs = fig.add_subplot(1,3,1, projection = wcs_proj)
    # defines coordinates of limits
    print "world ", world[0], world[1]
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
    print "pix ", pix[0], pix[1] 
    
    plt.show()

    return 0

if __name__ == '__main__':
    sys.exit(main())
