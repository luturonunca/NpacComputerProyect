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

def nearby_galaxies(fluxes, coords_input, coords_detected, wl):
    """
    A function that returns the number of input galaxies that are nearby to
    coord_detected.
    """

    for coord_in in (coords_input):
        for coord_de in (coords_detected):
            






#    nearby = []
#    n_nearby = []
#    nearby_coords = []
#    n_nearby_coords = []
#    r_ap = 2.5 * library.get_pixFWHM(wl)
#
#    for coord_in in (coords_input):
#        for coord_de in (coords_detected):
#            if abs(coord_de[0] - coord_in[0]) > r_ap\
#                                    or abs(coord_de[1] - coord_in[1]) > r_ap:
#                continue
#            else:
#                nearby.append(True)
#                nearby_coords.append(coord_in)
#    
#    for i in range(len(coords_input)):
#        for j in range(len(nearby_coords)):
#            if j == i:
#                continue
#            else:
#                n_nearby.append
#
#
#    if len(nearby) == len(n_nearby):
#        return nearby, nearby_coords
#    else:
#        #Merge arrays JOTA
#        for k in range(len(n_nearby))
#            #kkkkkkkk
#
#
#    return nearby, coords_input

#def nearby(coords_input, coords_detected, wl):
#    nearby = False
#    r_ap = 2.5 * library.get_pixFWHM(wl)
#    for coord_de in (coords_detected):
#        if abs(coord_de[0] - coords_input[0]) > r_ap\
#                                or abs(coord_de[1] - coords_input[1]) > r_ap:
#            continue
#        else:
#            nearby = True
#    return nearby

def main():
    fig = plt.figure()

    
    ### SHORT WAVELENGTH
    above_S, fits_data_S = library.get_data("S")
    image_S = above_S * fits_data_S
    _ , header_S = fits.getdata(get_filename("S"),header=True)
    # Create mask without NaN
    m_array_S = library.create_marray(fits_data_S)
   
    # check limits 289.345018791 -33.6099453587
    # Plotting section
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
    galaxy_counter = 0
    base = np.zeros((len(image_S[0]),len(image_S)))
    coord_array = []
    noise_me = 0
    clean_me = 0
    num_nans = 0
    diff = 0
    ## Arrays to plot photometric accuracy
    photo_accuracy = []
    flux_clean = []
    ## Initialize array for completeness plot it'll contain pairs of
    ## jansky, ratio.
    ratio = []
    flux_array = []
    matching_array = []

    while galaxy_counter < 1:
        counter = 0
        #################  Creates base arrays #######################      
        noise_base = np.copy(image_S)                                       
        base = np.zeros((len(image_S[0]),len(image_S))) 
        
        ## Local array of fluxes for the nine simulated galaxies below
        local_fluxes = []
        local_coords = []
        
        while counter < 9:
            counter += 1
            galaxy_counter += 1
            ################## Creates random Galaxies ####################
            #gauss_kernel = 5 * Gaussian2DKernel(2).array
            #witdth_rand = np.random.randint(1,3)
            gauss_kernel = np.random.uniform(0.1,20, size=1)[0] *\
                           Gaussian2DKernel(2).array
            coords_gal = np.random.randint(20,len(image_S[0])-20, size=2)
            coord_array.append(coords_gal) 
            local_coords.append(coords_gal) 
            # To plot completeness, we need flux value of sim galaxy
            local_fluxes.append(np.amax(gauss_kernel))
            flux_array.append(np.amax(gauss_kernel))
            
            for j in range(0, len(gauss_kernel)):
                for i in range(0, len(gauss_kernel[0])):
                    noise_base[coords_gal[1]+j-8][coords_gal[0]+i-8] += \
                               gauss_kernel[i][j]
            for j in range(0, len(gauss_kernel)):
                for i in range(0, len(gauss_kernel[0])):
                    base[coords_gal[1]+j-8][coords_gal[0]+i-8] += \
                                            gauss_kernel[i][j]
            ################## integrates both w noise n w'out ############
            integ_clean = library.photometrySimple(base,coords_gal,"S")
            integ_noise = library.photometrySimple(noise_base,coords_gal,"S")

            if math.isnan(integ_noise[1]) is False:
                noise_me += integ_noise[1]
                clean_me += integ_clean[1]
                diff += abs(integ_noise[1] - integ_clean[1])
            else:
                num_nans += 1
            
            ## Photometric accuracy part
            photo_accuracy.append((integ_clean[2] - integ_noise[2]) / integ_clean[2])
            flux_clean.append(integ_clean[2])
        
        
        ## Completeness part
        filt_n = library.filter_direct("S", noise_base)
        
		## Version of locator for completeness
        n_sig = 2
        sigma_n, mean_n = library.get_gauss_sigma(filt_n)
        mask = np.zeros((len(filt_n[0]), len(filt_n)))
        mask[filt_n >= n_sig * sigma_n + mean_n] =\
                 filt_n[filt_n >= n_sig * sigma_n + mean_n]

        suma_array = []
        centroides_array = []
        lumi_array = []                                                         
        #print "Length mask_plot ", len(mask),\
        #      "and mask_plot[1]", len(mask[1])                             
                                                                            
        for l in range(0, len(mask)):                                      
            for k in range(0, len(mask[1])):                               
                if mask[l][k] == 0:                                        
                    continue                                                    
                else:                                                           
                    mask, suma, ext, lumi = library.cluster_scan([k,l], mask) 
                    if suma < 5:                                                
                        continue                                                
                    else:                                                       
                        suma_array.append(suma)                                 
                        centroides_array.append(library.centroides(ext))       
                        lumi_array.append(lumi)
        
        
        matched, unmatched = nearby_galaxies(local_fluxes, local_coords, centroides_array, "S")
        print "List unmatched ", unmatched
        # little loop to append coherently in big array
        for m in range(0, len(matched)):
            matching_array.append(matched[m])
        print "Flux and matching: ", flux_array, matching_array

    print "Lengths of flux and matching: ", len(flux_array), len(matching_array)
    

#    print "\n## CENTROIDS ##\n Centroids and lumi :", centroides_array,\
#           lumi_array, "Length of centroids: ", len(centroides_array)
    #print filtered_noise
    #print base_noise
    
    #print "\n\n##Photo accuracy##\n\n Accuracy: ", photo_accuracy,\
    #      "Input flux: ", flux_clean
    
    print "Clean integral = ", clean_me /900., ",with noise = ", noise_me/900.
    print "Difference = ", diff/ 900.
    #print "pix ", pix[0], pix[1] 
    #ax_wcs.set_xlim(pix0[0], pix[0])
    #ax_wcs.set_ylim(pix0[1], pix[1])
    ax_wcs.imshow(noise_base, origin='lower', interpolation='None')
    ax_wcs1.imshow(base, origin='lower', interpolation='None')
    ax_wcs.coords['ra'].set_ticks(color='red')
    ax_wcs.coords['dec'].set_ticks(color='red')    
    plt.show()

    ## Photometric accuracy plots
    # Latex details
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    
    fig2 = plt.figure()
    ax = fig2.add_subplot(1,1,1)
    ax.scatter(flux_clean, photo_accuracy)
    ax.set_title('Photometric accuracy')
    plt.show()
    #plt.ylabel("$\frac{F_{input}-F_{measured}}{F_{input}}$")
    plt.ylabel("Normalized difference of flux")
    plt.xlabel("Flux")
    return 0

if __name__ == '__main__':
    sys.exit(main())
