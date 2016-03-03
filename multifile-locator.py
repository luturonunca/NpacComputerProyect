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
str_factor = (np.pi/180)**2
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
def compute_size_field(x,y,angle):
    size = x * y * (angle**2) * str_factor
    return size


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
    total_size_S = 0
    total_size_M = 0
    total_size_L = 0
    for j in range(len(files_array)):
        waverange =  files_array[j][1][-7]
        above, fits_data, header = get_data_path(waverange,files_array[j][0], header=True)
        image = above * fits_data
        total_size_S += compute_size_field(len(image[0]),len(image),\
                               abs(header['CDELT1']))
        print header['CDELT1']
        total_size_M += compute_size_field(len(image[0]),len(image),\
                               abs(header['CDELT1']))
        total_size_L += compute_size_field(len(image[0]),len(image),\
                               abs(header['CDELT1']))
        
        #print len(image),len(image[0]), abs(header['CDELT1'])
        dictionary = get_positions(files_array[j])
        coordinates = dictionary['centroids']
        for i in range(0, len(coordinates)):
            flux = photometrySimple(image,coordinates[i],waverange)
            if np.isnan(flux[2])== True:                                           
                 continue  
            if waverange == "S":
                cont += 1
                if np.isnan(flux[2])== True:
                    continue
                flux_array_S.append(flux[2])
            elif waverange == "M":
                cont += 1      
                if np.isnan(flux[2])== True:                                    
                    continue
                flux_array_M.append(flux[2])

            elif waverange == "L":                                                   
                cont += 1           
                if np.isnan(flux[2])== True:                                    
                    continue
                flux_array_L.append(flux[2])
    binboundaries = [0.02, 0.029,0.051,0.069,0.111,0.289,0.511] 
    bincenter = [0.0238, 0.0375, 0.0589, 0.0859, 0.1662, 0.3741]
    plt.close("all")
    plt.xscale("log")
    F_S,_ = np.histogram(flux_array_S, bins=binboundaries )
    print bincenter
    FF_S = []
    FF_M = []
    FF_L = []
    for i in range(len(F_S)):
            
            print float(F_S[i]), bincenter[i]**2.5, total_size_S
            FF_S.append(float( F_S[i])*(bincenter[i]**2.5)/float(total_size_S))  
    plt.plot(bincenter, FF_S, 'g^', markersize=20)
    
    F_M, _ = np.histogram(flux_array_M, bins=binboundaries )                     
    print bincenter                                                         
    for i in range(len(F_M)):
            print total_size_M
            FF_M.append(float( F_M[i])*(bincenter[i]**2.5)/(total_size_M))
    print FF_M        
    plt.plot(bincenter, FF_M, 'r^', markersize=20)


    F_L,_ = np.histogram(flux_array_L, bins=binboundaries )                     
    print bincenter                                                             
    for i in range(len(F_L)):     
            FF_L.append(float( F_L[i])*(bincenter[i]**2.5)/float(total_size_L)) 
    plt.plot(bincenter, FF_L, 'b^', markersize=20)
    
    plt.show()
    print FF_S
    print final_S, final_M, final_L

    return 0


if __name__ == '__main__':
    sys.exit(main())
