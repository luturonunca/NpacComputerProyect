#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import math
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import os
from library import *

plt.ion()
path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"         
filenameS= "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_PSW.fits" 
filenameM= "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_PMW.fits"
filenameL= "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_PLW.fits"
def main():

    mex_S = filter("S", path+filenameS)
    mex_M = filter("M", path+filenameM)
    mex_L = filter("L", path+filenameL)
    
    # Number of sigmas to be taken for the threshold
    n_sig = 2
    sigma_S, mean_S = get_gauss_sigma(mex_S)
    sigma_M, mean_M = get_gauss_sigma(mex_M)
    sigma_L, mean_L = get_gauss_sigma(mex_L)
    
    mask_S = np.zeros((len(mex_S[0]), len(mex_S)))
    mask_S[mex_S >= n_sig * sigma_S + mean_S] =\
                                    mex_S[mex_S >= n_sig * sigma_S + mean_S]

    mask_M = np.zeros((len(mex_M[0]), len(mex_M)))
    mask_M[mex_M >= n_sig * sigma_M + mean_M] =\
                                    mex_M[mex_M >= n_sig * sigma_M + mean_M]

    mask_L = np.zeros((len(mex_L[0]), len(mex_L)))
    mask_L[mex_L >= n_sig * sigma_L + mean_L] =\
                                    mex_L[mex_L >= n_sig * sigma_L + mean_L]

    pixels = None
    fig, main_axes = plt.subplots(1,3, figsize=(15,5))
    main_axes[0].imshow(mask_S, origin="lower", interpolation="None")
    main_axes[1].imshow(mask_M, origin="lower", interpolation="None")
    main_axes[2].imshow(mask_L, origin="lower", interpolation="None")
    plt.show()


    d_gal = {}

    for mask in [(mask_S,"S"), (mask_M,"M"), (mask_L,"L")]:
        suma_array = []
        centroides_array = []
        lumi_array = []
        mask_plot = np.copy(mask[0])
        print "Length mask_plot ", len(mask_plot),\
              "and mask_plot[1]", len(mask_plot[1])
    
        for j in range(0, len(mask_plot)):
            for i in range(0, len(mask_plot[1])):
                if mask_plot[j][i] == 0:
                    continue
                else:
                    #print "Coordinates ", i, j, "Value ", mask_plot[j][i]
                    mask_plot, suma, ext, lumi = cluster_scan([i,j], mask_plot)
                    if suma < 5:
                        continue
                    else:
                        suma_array.append(suma)
                        centroides_array.append(centroides(ext))
                        lumi_array.append(lumi)
        
        
        d_gal[mask[1]] = {"centroids": centroides_array, "n_pix": suma_array,\
                          "luminosity_first": lumi_array}

        #print "Suma\n", suma_array, "\nLumi\n", lumi_array, "\nCentroids\n",\
        #       centroides_array
        #print "\nDICTIONARY\n", d_gal["S"]["centroids"][1]

        # Safety check
        #print "Lengths: ", len(suma_array), len(lumi_array), len(centroides_array)

    return 0


if __name__ == '__main__':
    sys.exit(main())
