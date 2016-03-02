#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test and simulation to validate photometry meaassurement

"""
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

def photometryT(data,centroid,wl):
    """
    Test version of photometry function
    recieves centroid in pix cordinates
    makes de integral of the galaxy centered at
    centroid
    """
    #wcs = WCS(header)
    fwhm = library.get_pixFWHM(wl)
    #cent_pix = wcs.wcs_world2pix(centroid[0], centroid[1],0)
    cent_pix = centroid
    r_ap = 2.5 * fwhm
    r_int = 3.5 * fwhm
    r_ext = 4.5 * fwhm
    # set iteration limits
    # x_l -> leftend; x_r -> rightend
    # y_u -> upend; y_b -> bottomend
    x_l = int(cent_pix[0] - r_ext)
    x_r = int(cent_pix[0] + r_ext)
    y_u = int(cent_pix[1] + r_ext)
    y_b = int(cent_pix[1] - r_ext)
    integral = 0
    background = 0
    size_ring = 0
    size_bump = 0
    log_bag = 0
    for j in range(y_b, y_u):
        for i in range(x_l, x_r):
            x_local = i - cent_pix[0]
            y_local = j - cent_pix[1]
            r = np.sqrt((x_local**2)+(y_local**2))

            # Check that values are inside array
            if i >= len(data[0]) or i < 0:
                continue
            if j >=len(data) or j < 0:
                continue
            if r <= r_ap:
                size_bump += 1
                integral +=   data[j][i]
            elif r >= r_int and r <= r_ext:
                size_ring += 1.
                background += data[j][i]
                log_bag += math.pow(10,(data[j][i]/10.)) 
    norm_bg = background / size_ring
    log_ave = np.log10(log_bag / size_ring)
    reduced_sig = integral - (size_bump * norm_bg)
    return  [integral, reduced_sig, size_bump], [(size_bump * norm_bg) , norm_bg, log_ave]



def main():
    # firs test    
    array = np.zeros((500,500)) 
    array[250][250] = 10
    result = photometryT(array,[250,250],"S")
    print "## First test zeros array with one 10 in de center ##"
    print "integrals: circle = ", result[0][0], ", normalized circle = ", result[0][1]
    print "backgrouns: in cirque = ", result[1][0], ", average per pixel = ", result[1][1]

    
    # Second test
    array2 = np.zeros((500,500))
    array2.fill(100)
    array2[250][250] += 10
    result2 = photometryT(array2,[250,250],"S")
    print "\n"
    print "## Second test constant array  with one constant+10 in de center ##"
    print "integrals: circle = ", result2[0][0], ", normalized circle = ", result2[0][1]
    print "backgrouns: in circle = ", result2[1][0], ", average per pixel = ", result2[1][1]


    # Third test
    array3 = np.random.uniform(0,3, size=(300, 300))
    array3[150][150] += 10
    result3 = photometryT(array3,[150,150],"S")

    print "\n"
    print "## thrid test random array beetwen 0 and 3 with one 10 in de center ##"
    print "integrals:     circle = ", result3[0][0], ", normalized circle = ", result3[0][1]
    print "backgrouns: in circle = ", result3[1][0], ", average per pixel = ", result3[1][1]
    final_log = 0 
    final = 0
    cont = 0
    while cont <10000.:
        array3 = np.random.uniform(0,3, size=(300, 300))
        array3[150][150] += 10
        result3 = photometryT(array3,[150,150],"S")
        final += result3[0][1]
        final_log += result3[0][0]-(result3[1][2] * result3[0][2])
        cont +=1.
    print "\n#################### After 10,000 calculations ####################  "
    print "average in circle integral = ", final/cont
    
    return 0

if __name__ == '__main__':
    sys.exit(main())

