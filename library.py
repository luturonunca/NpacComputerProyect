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
    """
    Gets Full Width Maximun for the diferent images in arcseconds
    """
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
    Gets Full Width Maximun for the diferent images in pixels
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

def filter(wavelength_range,filepath):                                                   
    #print wavelength_range                                                     
    above, fits_data = get_data_path(wavelength_range,filepath)                       
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
    fit, covariant = curve_fit(gaussian_fit, normal_x, normal_y, maxfev=2000)        
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
    if x < (len(data[y])-1) and data[y][x+1] > 0:                               
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

def photometry(data,centroid,wl,header):
    """
    recieves centroid in ra dec
    makes de integral of the galaxy centered at 
    centroid
    """
    wcs = WCS(header)
    fwhm = get_pixFWHM(wl)
    cent_pix = wcs.wcs_world2pix(centroid[0], centroid[1],0)
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
    norm_bg = background / size_ring
    reduced_sig = integral - (size_bump * norm_bg)
    return  [integral, reduced_sig]

def photometrySimple(data,centroid,wl):
    """
    Test version of photometry function
    recieves centroid in pix cordinates
    makes de integral of the galaxy centered at
    centroid
    """
    #wcs = WCS(header)
    fwhm = get_pixFWHM(wl)
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
    array_circle = []
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
                array_circle.append(data[j][i])
                size_bump += 1
                integral +=   data[j][i]
            elif r >= r_int and r <= r_ext:
                size_ring += 1.
                background += data[j][i]
                log_bag += math.pow(10,(data[j][i]/10.))
    norm_bg = background / size_ring
    log_ave = np.log10(log_bag / size_ring)
    reduced_sig = integral - (size_bump * norm_bg)
    flux_J = np.amax(array_circle)- norm_bg
    return  [integral, reduced_sig,flux_J]


def get_data_path(W,filepath, header = False):                                                                
      """                                                                         
      gets either S, M or L as string representing the dataset wavelengh          
      returns: the cutted data, i.e centered, and a boolean matrix                
      that represent a threathole in the coverage matrix.                         
                                                                                 
      """                                                                         
                                                                                  
      # Load the data                                                             
      fits_data_uncut, fits_header = fits.getdata(filepath,"image", header=True)                       
                                                                                  
      fits_coverage_uncut = fits.getdata(filepath,"coverage")                
      # Define cuts                                                               
      ycenter = int(len(fits_data_uncut)/2)                                       
      xcenter = int(len(fits_data_uncut[0])/2)                                    
      cut = int(ycenter/2)                                                        
      # Cut array to focus on center                                              
      fits_data = fits_data_uncut[xcenter-cut:xcenter+cut,ycenter-cut:ycenter+cut]
      fits_coverage = fits_coverage_uncut[xcenter-cut:xcenter+cut,ycenter-cut:ycenter+cut    ]
                                                                                  
      # Create masked array from fits_data                                        
      m_array = create_marray(fits_data)                                          
      # Create masked array from fits_data maskin with coverage                   
      mask_coverage = m_array.mask * fits_coverage                                
      above = fits_coverage > 8                                                   
      ## build this outside 
      # image = (above * fits_data)                                                 
      if header == True:
        return above, fits_data, fits_header
      else:
        return above, fits_data


