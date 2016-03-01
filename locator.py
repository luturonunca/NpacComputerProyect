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

def main():

    mex_S = filter("S")
    mex_M = filter("M")
    mex_L = filter("L")
    
    sigma_S, mean_S = get_gauss_sigma(mex_S)
    sigma_M, mean_M = get_gauss_sigma(mex_M)
    sigma_L, mean_L = get_gauss_sigma(mex_L)
    
    mask_S = np.zeros((len(mex_S[0]), len(mex_S)))
    mask_S[mex_S >= 3 * sigma_S + mean_S] = mex_S[mex_S >= 3 * sigma_S + mean_S]

    mask_M = np.zeros((len(mex_M[0]), len(mex_M)))
    mask_M[mex_M >= 3 * sigma_M + mean_M] = mex_M[mex_M >= 3 * sigma_M + mean_M]

    mask_L = np.zeros((len(mex_L[0]), len(mex_L)))
    mask_L[mex_L >= 3 * sigma_L + mean_L] = mex_L[mex_L >= 3 * sigma_L + mean_L]

    pixels = None
    fig, main_axes = plt.subplots(1,3, figsize=(15,5))
    main_axes[0].imshow(mask_S, origin="lower", interpolation="None")
    main_axes[1].imshow(mask_M, origin="lower", interpolation="None")
    main_axes[2].imshow(mask_L, origin="lower", interpolation="None")
    plt.show()

    
    
    suma_array = []
    centroides_array = []
    lumi_array = []
    mask_plot = np.copy(mask_S)
    print "Length mask_plot ", len(mask_plot), "and mask_plot[0]", len(mask_plot[0])

    for j in range(0, len(mask_plot) - 1):
        for i in range(0, len(mask_plot[1]) - 1):
            if mask_plot[j][i] > 0:
                #print "Coordinates ", i, j, "Value ", mask_plot[j][i]
                mask_plot, suma, extrema, lumi = cluster_scan([i,j], mask_plot)
                suma_array.append(suma)
                centroides_array.append(centroides(extrema))
                lumi_array.append(lumi)
            else:
                continue

    print "Suma\n", suma_array, "\nLumi\n", lumi_array, "\nCentroids\n",\
           centroides_array
    print "Lengths: ", len(suma_array), len(lumi_array), len(centroides_array)

    return 0


if __name__ == '__main__':
    sys.exit(main())
