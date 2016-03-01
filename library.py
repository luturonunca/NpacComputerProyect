#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
personal library
"""
import sys                                                                      
import math                                                                     
import numpy as np                                                              
import numpy.ma as ma                                                           
import matplotlib.pyplot as plt                                                 
from astropy.io import fits
from astropy.wcs import WCS                                                     
from astropy.convolution import Gaussian2DKernel, MexicanHat2DKernel, convolve
import os
from scipy.optimize import curve_fit

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

def filter(wavelength_range):                                                   
    #print wavelength_range                                                     
    above, fits_data = get_data(wavelength_range)                       
    image = above * fits_data                                                   
                                                                                
    # Create mask without NaN                                                   
    m_array = create_marray(fits_data)                                  
                                                                                
    # Make the Mexican hat kernel and convolve                                  
    # The values given to the kernels should be 1.7328, 2.3963, 3.5270          
    # for S, M and L.                                                           
    mex = MexicanHat2DKernel(get_pixFWHM(wavelength_range)) 
    mex_convol = convolve(image, mex, boundary = 'extend')                      
    m_mex = ma.masked_array(mex_convol, ~m_array.mask)                          
                                                                                
    # Make the gaussian kernel and convolve                                     
    gauss = Gaussian2DKernel(stddev=get_pixFWHM(wavelength_range)) 
    gauss_convol = convolve(image, gauss, boundary='extend')                    
    m_gauss= ma.masked_array(gauss_convol, ~m_array.mask)                       
    #c_gauss = np.multiply(z, m_array.mask)                                     
                                                                                
#    ### Uncomment the next few lines to produce plots for the filtered mexican 
#    ### hat, gaussian and raw data. 151,173s/^/#/g                             
#                                                                               
#                                                                               
#    # Plot the figures; the min and the max values are reset in some of them   
#    # for visualization purposes.                                              
#    pixels = None                                                              
#    fig, main_axes = plt.subplots(1,3, sharex=True, sharey=True, figsize=(15,5))
#    #main_axes[0][0].imshow(mex, origin="lower", interpolation="None")         
#    #main_axes[0][0].set_title("MexHat kernel")                                
#                                                                               
#    main_axes[0].imshow(m_mex,\                                                
#                           origin="lower", interpolation="None")               
#    main_axes[0].set_title("Convolution with MexHat")                          
#                                                                               
#    main_axes[1].imshow(m_gauss, origin="lower", interpolation="None")         
#    main_axes[1].set_title("Convolution with gaussian")                        
#                                                                               
#                                                                               
#    main_axes[2].imshow(image, origin="lower", interpolation="None")           
#    main_axes[2].set_title("Data")                                             
#                                                                               
#    plt.show()                                                                 
                                                                                
#    ### Uncomment the next few lines to produce and print statistical values   
#    mean_data = (np.nanmax(fits_data) - np.nanmin(fits_data))/2                
#    print "Mean", mean_data, "\nMax", np.nanmax(fits_data),\                   
#          "\nMin", np.nanmin(fits_data)
    return m_mex


def gaussian_fit(w, max_value, mean_value, sigma):
    """
    Gaussian fit
    """
    # pylint: disable=E1101
    z = max_value*np.exp(-((w-mean_value)**2)/(2*(sigma**2)))
    return z

def get_gauss_sigma(data):                                                      
    """                                                                         
    Modified function from library. Does the same but returns more              
    data from fit                                                               
    """                                                                         
    bin_values, bin_boundaries = np.histogram(data, 200)                        
    my = np.float(np.max(np.nonzero(bin_values)))                               
    mx = np.float(np.max(bin_boundaries))                                       
    normal_y = bin_values/my                                                    
    normal_x = bin_boundaries[:-1]/mx                                           
    x = normal_x * mx                                                           
    y = normal_y * my                                                           
    fit, covariant = curve_fit(gaussian_fit, normal_x, normal_y)        
    maxvalue = fit[0] * my                                                      
    background = fit[1] * mx                                                    
    dispersion = fit[2] * mx                                                    
    return abs(dispersion), background

def check_n_move(coords, data, suma, lumi):                                     
    """                                                                         
    check neighbours                                                            
    """                                                                         
    x = coords[0]                                                               
    y = coords[1]                                                               
    # check right                                                               
    if data[y][x+1] > 0 and x < (len(data[y])-1):                               
        # print "moved right"                                                   
        suma += 1                                                               
        lumi += data[y][x+1]                                                    
        data[y][x+1] = 0                                                        
        x += 1                                                                  
    # check top                                                                 
    elif  y < (len(data)-1) and data[y+1][x] > 0:                               
        # print "moved down"                                                    
        suma += 1                                                               
        lumi += data[y+1][x]                                                    
        data[y+1][x] = 0                                                        
        y += 1                                                                  
    # check left                                                                
    elif data[y][x-1] > 0 and x != 0:                                           
        # print "moved left"                                                    
        suma += 1                                                               
        lumi += data[y][x-1]                                                    
        data[y][x-1] = 0                                                        
        x = x-1                                                                 
    # check down                                                                
    elif data[y-1][x] > 0 and x != 0:                                           
        # print "moved up"                                                      
        suma += 1                                                               
        lumi += data[y-1][x]                                                    
        data[y-1][x] = 0                                                        
        y = y-1                                                                 
    else:                                                                       
        # print "stayed                                                         
        x = coords[0]                                                           
        y = coords[1]                                                           
    return x, y, data, suma, lumi

def cluster_scan(initial, dat):                                                 
    """                                                                         
     scan cluster until  finish                                                 
    """                                                                         
    # Extremes of the cluster [x_min, x_max, y_min(fixed), y_max]
    extremos = [initial[0], initial[0], initial[1], initial[1]]                 
    dat[initial[1]][initial[0]] = False                                         
    camino = [initial]                                                          
    current1 = initial                                                          
    suma = 1                                                                    
    lumi = 0                                                                    
    complete = False                                                            
    while complete is False:                                                    
        x, y, dat, suma, lumi = check_n_move(current1, dat, suma, lumi)         
        if x > extremos[1]:                                                     
            extremos[1] = x                                                     
        elif x < extremos[0]:                                                   
            extremos[0] = x                                                     
        elif y > extremos[3]:                                                   
            extremos[3] = y                                                     
                                                                                
        if [x, y] == camino[-1]:                                                
            if suma == 1:                                                       
                # -flag- print "unborn"                                         
                complete = True                                                 
            elif len(camino) > 1:                                               
                # -flag- print "sigue, camino:", len(camino)                    
                [x, y] = camino[-2]                                             
                del camino[-1]                                                  
                # -flag- print "sigue, camino:", len(camino)                    
                                                                                
            elif len(camino) == 1:                                              
                # -flag- print "finish"                                         
                complete = True                                                 
        else:                                                                   
            # -flag- print "appending "                                         
            camino.append([x, y])                                               
        current1 = [x, y]                                                       
        #if suma==100:                                                          
        #    complete = True                                                    
    return dat, suma, extremos, lumi

def centroides(arreglo):                                                        
    x_centroid = ((arreglo[1] - arreglo[0])/2.+ arreglo[0])                     
    y_centroid = ((arreglo[3] - arreglo[2])/2. + arreglo[2])                    
    return x_centroid, y_centroid


def create_marray(image):
    """
    gets data to mask
    creates mask without NaN
    returns mask
    """
    mask = np.zeros((len(image), len(image[0])), dtype=bool)
    for j in range (len(image)):
        for i in range (len(image[j])):
            if math.isnan(image[j][i]):
                continue
            else:
                mask[j][i] = True

    # Create the masked array. 1-mask, be careful.
    m_array = ma.masked_array(image, mask)
    return m_array


def get_data(W):
    """
    gets either S, M or L as string representing the dataset wavelengh
    returns: the cutted data, i.e centered, and a boolean matrix
    that represent a threathole in the coverage matrix.
    
    """
    
    path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
    filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_P"\
    +W+"W.fits"
    # Load the data
    fits_data_uncut = fits.getdata(path+filename,"image")
    
    fits_coverage_uncut = fits.getdata(path+filename,"coverage")
    # Define cuts
    ycenter = int(len(fits_data_uncut)/2)
    xcenter = int(len(fits_data_uncut[0])/2)
    cut = int(ycenter/2)
    print "centers ", float(xcenter)/ycenter
    # Cut array to focus on center
    fits_data = fits_data_uncut[xcenter-cut:xcenter+cut,ycenter-cut:ycenter+cut]
    fits_coverage = fits_coverage_uncut[xcenter-cut:xcenter+cut,ycenter-cut:ycenter+cut]

    # Create masked array from fits_data
    m_array = create_marray(fits_data)
    # Create masked array from fits_data maskin with coverage
    mask_coverage = m_array.mask * fits_coverage
    above = fits_coverage > 8

    image = (above * fits_data)
    
    return above, fits_data



def transform(array,W):
    path = "/home/abeelen/Herschel/DDT_mustdo_5/PLCK_SZ_G004.5-19.6-1/"
    filename = "OD1271_0x50012da8L_SpirePhotoLargeScan_PLCK_SZ_G004.5-19.6-1_destriped_P"\
    +W+"W.fits"
    header = fits.getheader(path+filename)
    w = WCS(header)
    array2 = w.wcs_pix2world(array,0)
    return array2


