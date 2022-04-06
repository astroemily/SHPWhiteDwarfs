#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 12:53:09 2022

@author: emilyroberts
"""

import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd

line = []
line_no = 0
name = [] #empty list for star ID number
colour = [] #empty list for star colour R-B
abs_mag_b = [] #empty list for absolute magnitudes in B band
distance = [] #empty list for distance to star
ruwe_data = []
excess_data = []

wd_only_mag = []
wd_only_colour = []
wd_only_names = []
    
  
table = []
 
stars = 0
with open('gaia_edr3_250pc.txt', 'r') as data_file:

    for line in data_file:
        data = line.split() #split line of file where there is a space
        stars +=1
        try:
            ruwe = float(data[14]) #some of the data lines don't seem to have these two statistics on the line
        except:  
            #need to skip those lines!
            print("skipped a line 14")
            continue
        else:
            #print("ruwe ", ruwe)
            ruwe = float(data[14])
            phot_excess = float(data[15])
            ruwe_data.append(ruwe)
            excess_data.append(phot_excess)
        
        
        
            star_name = data[0]
        
            star_blue = data[12]
            star_red = data[13] #separating out the elements using the read me nick gave for which values needed
            parallax = data[5]
        
            star_blue_int = float(star_blue) #turning them back into integers not strings
            star_red_int = float(star_red)
            parallax_int = float(parallax)
            
            star_colour = star_blue_int - star_red_int #calculating star colour
            #print(star_name, "has colour ", colour)
            colour.append(star_colour)
            
            #got the colour, now need the absolute magnitude in b band for HR diagram
            parallax_int_pc = parallax_int/1000 #converting parallax to arcsecs from milli arc secs
            distance_pc = 1/parallax_int_pc #distance in parsecs to star
        
    
            abs_mag = star_blue_int - 5*np.log10(distance_pc/10) #distance  modulus calculation for absolute mag
            #print("absolute b band magnitude = ", abs_mag)
"""            
plt.plot(colour, excess_data, 'ro', markersize=1, linestyle='None')
plt.axhline(1.1,color='black')
plt.title("Photo Excess Value against Colour for Full Gaia Data")
plt.xlabel("Colour")
plt.ylabel("Photo Excess Value")
plt.ylim(0,5)

plt.show()


plt.plot(colour, ruwe_data, 'bo', markersize=1, linestyle='None')
plt.axhline(1.1, color='black')
plt.title("Renormalised Unit Weight Error against Colour for Full Gaia Data")
plt.xlabel("Colour")
plt.ylabel("RUWE")
plt.ylim(0,5)

plt.show()
"""
x = colour
y = ruwe_data
z = excess_data


x_min = np.min(x)
x_max = np.max(x)

y_min = np.min(y)
y_max = np.max(y)

x_bins = np.linspace(x_min, x_max, 10000)
y_bins = np.linspace(y_min, y_max, 10000)

plt.hist2d(x, y, bins =[x_bins, y_bins],norm=mcolors.PowerNorm(0.3),cmap='inferno')
plt.title('BP-RP colour against RUWE')
plt.axhline(1.1, color='white')
plt.ylim(np.min(y), 5)
plt.xlabel("BP-RP colour")
plt.ylabel("RUWE value")
plt.xlim(-1,5)
plt.colorbar()



z_min = np.min(z)
z_max = np.max(z)

# z_bins = np.linspace(z_min,z_max, 10000)
# z_line15 = [1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1]
# z_line2 = [1,1.2,1.4,1.6,1.8,2,2.2,2.4]
# z_line216 = [0.984,1.2,1.416,1.632,1.848,2.064,2.28,2.496]
# colour_line = [-1,0,1,2,3,4,5,6]
# plt.hist2d(x, z, bins =[x_bins, z_bins],norm=mcolors.PowerNorm(0.3),cmap='inferno')
# plt.title('BP-RP colour against Photo Excess')
# plt.plot(colour_line,z_line15, color='white')
# plt.ylim(np.min(z), 3)
# plt.xlabel("BP-RP colour")
# plt.ylabel("Photo Excess value")
# plt.xlim(-1,5)
# plt.colorbar()
                     
                     
                     