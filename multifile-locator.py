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
def file_traveler(path):
    """
    iterates over directory on two levels and
    returns an array of full path+filename
    """
    directory = os.listdir(path)
    files_array = []
    for fn in directory:
        cop = os.listdir(path+fn)                                                   
        for files in cop:                                                           
            files_array.append([path+fn+"/"+files,files])
    return files_array



def get_positions(fileinfo):
    mex = filter(fileinfo[1][-7],fileinfo[0])
    
    # Number of sigmas to be taken for the threshold
    n_sig = 3
    sigma, mean = get_gauss_sigma(mex)
    
    mask= np.zeros((len(mex[0]), len(mex)))
    mask[mex >= n_sig * sigma + mean] =\
                                    mex[mex >= n_sig * sigma + mean]


    pixels = None
    #fig, main_axes = plt.subplots(1,3, figsize=(15,5))
    #main_axes[0].imshow(mask_S, origin="lower", interpolation="None")
    plt.show()


    d_gal = {}
    suma_array = []
    centroides_array = []
    lumi_array = []
    mask_plot = np.copy(mask)
    #print "Length mask_plot ", len(mask_plot),\
              #"and mask_plot[1]", len(mask_plot[1])
    
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
        
        
    d_gal = {"centroids": centroides_array, "n_pix": suma_array,\
                      "luminosity_first": lumi_array}
    return d_gal

def main():
    ########## define directory where files are ####################
    path = "/home/abeelen/Herschel/DDT_mustdo_5/"                               
    files_array = file_traveler(path) #give [full-path-to-file,filename]
    flux_array_S = []
    flux_array_M = []
    flux_array_L = []
    x_array = []
    cont = 0
    final_S = [0,0,0,0,0,0]
    final_M = [0,0,0,0,0,0]
    final_L = [0,0,0,0,0,0]
    for i in range(len(files_array)):
        waverange =  files_array[i][1][-7]
        above, fits_data = get_data_path(waverange,files_array[i][0])
        image = above * fits_data
        dictionary = get_positions(files_array[i])
        coordinates = dictionary['centroids']

        for i in range(len(coordinates)):
            flux = photometrySimple(image,coordinates[i],waverange)
            if waverange == "S":
                cont += 1
                flux_array.append(flux[2])
                x_array.append(cont) 
                if flux[2] > 0.02 and flux[2] < 0.029:
                    final[0]+=1
                if flux[2] > 0.029 and flux[2] < 0.051:
                    final[1]+=1
                if flux[2] > 0.051 and flux[2] < 0.69:                           
                    final[2]+=1                                                 
                if flux[2] > 0.69 and flux[2] < 0.111:                           
                    final[3]+=1
                if flux[2] > 0.111 and flux[2] < 0.289:                          
                    final[4]+=1                                                
                if flux[2] > 0.289 and flux[2] < 0.511:                          
                    final[5]+=1

    return 0


if __name__ == '__main__':
    sys.exit(main())
