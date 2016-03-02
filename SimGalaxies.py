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

    
    ### SHORT WAVELENGH
    above_S, fits_data_S = library.get_data("S")
    image_S = above_S * fits_data_S
    _ , header_S = fits.getdata(get_filename("S"),header=True)
    # Create mask without NaN
    m_array_S = library.create_marray(fits_data_S)
   
    # check limits 289.345018791 -33.6099453587
    # Ploting section
    wcs_proj = WCS(header_S)
    world = wcs_proj.wcs_pix2world(57.,57., 0)                                
    zero = wcs_proj.wcs_pix2world(0.,0., 0)
    ax_wcs = fig.add_subplot(1,2,1, projection = wcs_proj)
    ax_wcs1 = fig.add_subplot(1,2,2, projection = wcs_proj)
    # defines coordinates of limitis
    # -flag- print "world ", world[0], world[1]
    pix = wcs_proj.wcs_world2pix(world[0],world[1], 0)
    pix0 = wcs_proj.wcs_world2pix(zero[0],zero[1], 0)
    
    ## simulating galaxies
    cont = 0
    base = np.zeros((len(image_S[0]),len(image_S)))
    while cont < 9:
        cont += 1
        gauss_kernel = np.random.uniform(1,3, size=1)[0] * Gaussian2DKernel(2).array
        coords_gal = np.random.randint(20,len(image_S[0])-20, size=2)
 
        #coords_gal = int(len(image_S)/2) , int(len(image_S)/2)
        for j in range(0, len(gauss_kernel)):
            for i in range(0, len(gauss_kernel[0])):
                image_S[coords_gal[1]+j-8][coords_gal[0]+i-8] += gauss_kernel[i][j]
        for j in range(0, len(gauss_kernel)):
                for i in range(0, len(gauss_kernel[0])):
                    base[coords_gal[1]+j-8][coords_gal[0]+i-8] += gauss_kernel[i][j]



    integ_clean = library.photometrySimple(base,[51,51],"S")
    coordss = [coords_gal[0]-8,coords_gal[1]-8]
    integ_noise = library.photometrySimple(image_S,coords_gal,"S")
    print "perfect integral = ", integ_clean[1], ",with noise = ", integ_noise[1]

    #print "pix ", pix[0], pix[1] 
    #ax_wcs.set_xlim(pix0[0], pix[0])
    #ax_wcs.set_ylim(pix0[1], pix[1])
    ax_wcs.imshow(image_S, origin='lower', interpolation='None')
    ax_wcs1.imshow(base, origin='lower', interpolation='None')
    ax_wcs.coords['ra'].set_ticks(color='red')
    ax_wcs.coords['dec'].set_ticks(color='red')    
    plt.show()

    return 0

if __name__ == '__main__':
    sys.exit(main())
