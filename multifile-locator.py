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
    n_sig = 2
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
                if suma < 1:
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
    
    number_of_Sfields = 0
    number_of_Mfields = 0
    number_of_Lfields = 0
    for j in range(len(files_array)):
        # check range of view of file
        waverange =  files_array[j][1][-7]
        # gets data and header
        above, fits_data, header = get_data_path(waverange,files_array[j][0],\
                                   header=True)
        print header['CDELT1'], files_array[j][1][-7]
        image = above * fits_data
        ### gets de dictionary and the coordinates of galaxies ###
        ## fot the current file but first ignore problematic file
        dictionary = get_positions(files_array[j])                              
        coordinates = dictionary['centroids']
        ### SHORT WAVELENGTH ###
        
        if waverange == "S":
            number_of_Sfields += 1
            ###### sums the size of the field that is observed in sr #####
            total_size_S += compute_size_field(len(image[0]),len(image),\
                               abs(header['CDELT1']))
            print compute_size_field(len(image[0]),len(image),abs(header['CDELT1']))
            for i in range(0, len(coordinates)):
                flux = photometrySimple(image,coordinates[i],waverange)
                if np.isnan(flux[2])== True:                                    
                    continue                                                    
                flux_array_S.append(flux[2])
        
        ### MEDIUM WAVELENGTH ### 
        elif waverange == "M":
            print "coco", number_of_Mfields
            number_of_Mfields += 1
            total_size_M += compute_size_field(len(image[0]),len(image),\
                               abs(header['CDELT1']))
            for i in range(0, len(coordinates)):                                
                flux = photometrySimple(image,coordinates[i],waverange)
                if np.isnan(flux[2])== True:                                    
                    continue                                                    
                flux_array_M.append(flux[2])
        
        ### LONG WAVELENGTH ### 
        elif waverange == "L":
            number_of_Lfields += 1
            total_size_L += compute_size_field(len(image[0]),len(image),\
                               abs(header['CDELT1']))
            print "imagened",len(image[0]),len(image)
            for i in range(0, len(coordinates)):                                
                flux = photometrySimple(image,coordinates[i],waverange)
                if np.isnan(flux[2])== True:                                    
                    continue                                                    
                flux_array_L.append(flux[2])
    
    ## literarute limits for histogram   
    binboundaries = [0.02, 0.029,0.051,0.069,0.111,0.289,0.511] 
    bincenter = [0.0238, 0.0375, 0.0589, 0.0859, 0.1662, 0.3741]
    ## completenes correction
    completeness = [0.5207333, 0.91823224, 0.98095451, 0.98438088,\
                    0.98823266, 0.99232246]
    completenessM = [0.50523013, 0.825, 0.925, 0.94548287,\
                    0.96344454, 0.97344334]
    completenessL = [0.514862, 0.60981912, 0.71747212, 0.83818626,\
                    0.87421197, 0.884297]
    ## plot secction
    #plt.close("all")
    fig = plt.figure()
    plt.ylabel(r'$ N^3$',fontsize=16)    
    plt.xlabel("Flux Density [Jy]")
    plt.axis('off')
    F_S,_ = np.histogram(flux_array_S, bins=binboundaries )
    print F_S
    FF_S = []
    FF_M = []
    FF_L = []
    ####### SHORT WAVELENGTH #######
    yerr1 = np.sqrt(F_S)
    for i in range(len(F_S)):
            yerr1[i]=(float(yerr1[i])*(bincenter[i]**2.5)/\
                        (float(total_size_S)))
            FF_S.append(float( F_S[i])*(bincenter[i]**2.5)*completeness[i]/\
                       (float(total_size_S)))  
    ax = fig.add_subplot(3,1,1)
    ## note that this is a relative error 
    compleatness_err = [65.874457619525543,87.475341287321271, \
                        90.413610240199276,90.571375168979301, \
                        90.748400353211011,90.935987559015004]
    yerr= [0,0,0,0,0,0]
    for m in range(len(completeness)):
            yerr[m] = np.sqrt(yerr1[m]**2 + compleatness_err[m]**2)
    
    ax.errorbar(bincenter[:-1], FF_S[:-1], yerr=yerr[:-1], fmt='ko', markersize=7)
    #ax.plot(bincenter, FF_S, 'ko', markersize=10)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.ylabel(r'$dN/dSd\Omega x S^{2.5}$[sr$^{-1}$Jy$^{1.5}$]')
    ax.text(0.3, 1600, r'250$\mu$')
    #ax.set_xticks([])
   


    ####### MEDIUM WAVELENGTH ####### 
    F_M, _ = np.histogram(flux_array_M, bins=binboundaries )                     
    yerr2 = np.sqrt(F_M)
    
    print bincenter     

    for i in range(len(F_M)):
            yerr2[i]=(float(yerr2[i])*(bincenter[i]**2.5)/\
                        (float(total_size_M)))
            FF_M.append(float( F_M[i])*(bincenter[i]**2.5)*completeness[i]/\
                        (total_size_M))
    ax2 = fig.add_subplot(3,1,2) 
    completenessM_err =[45.881647111526703, 58.630196997792872,\
                        62.081935107297248, 62.765531610377785,\
                        63.358916104996617, 63.686842571026133]
    ax2.plot(bincenter, FF_M, 'ko', markersize=10)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    plt.ylabel(r'$dN/dSd\Omega x S^{2.5}$[sr$^{-1}$Jy$^{1.5}$]')
    ax2.text(0.3, 150, r'350$\mu$')
    ax.set_xticklabels([])
    ax2.set_xticklabels([])
    xticklabels = ax.get_xticklabels() + ax2.get_xticklabels()
    plt.setp(xticklabels, visible=False)
    


    F_L,_ = np.histogram(flux_array_L, bins=binboundaries )                     
    print bincenter                                                             
    for i in range(len(F_L)):     
            FF_L.append(float( F_L[i])*(bincenter[i]**2.5)/float(total_size_L)) 
    ax3 = fig.add_subplot(3,1,3)  
    ax3.plot(bincenter, FF_L, 'ko', markersize=10)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.text(0.3, 1500, r'500$\mu$')
    plt.ylabel(r'$dN/dSd\Omega x S^{2.5}$[sr$^{-1}$Jy$^{1.5}$]')
    plt.xlabel('Flux Density [Jy]')
    plt.subplots_adjust(hspace = 0,
                        wspace = 0)
    plt.show()
    print FF_S
    print final_S, final_M, final_L

    return 0


if __name__ == '__main__':
    sys.exit(main())
