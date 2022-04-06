#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 14:51:32 2022

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
colour = [] #empty list for star colour
abs_mag_g = [] #empty list for absolute magnitudes in B band
distance = [] #empty list for distance to star
ruwe_data = []
excess_data = []

wd_only_mag = [] #lists for when i've filtered to keep only white dwarfs
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
            ruwe = float(data[14])
            phot_excess = float(data[15]) #if the two values exist, append them to lists
            ruwe_data.append(ruwe)
            excess_data.append(phot_excess)
        
        
        
            star_name = data[0]
        
            star_blue = data[12]
            star_g = data[11]
            star_red = data[13] #separating out the elements needed from the file
            parallax = data[5]
        
            star_blue_int = float(star_blue) #turning them back into integers not strings
            star_red_int = float(star_red)
            parallax_int = float(parallax)
            star_g_int = float(star_g)
            
            star_colour = star_blue_int - star_red_int #calculating star colour
            
            #got the colour, now need the absolute magnitude in g band for HR diagram
            
            parallax_int_pc = parallax_int/1000 #converting parallax to arcsecs from milli arc secs
            distance_pc = 1/parallax_int_pc #distance in parsecs to star
        
    
            abs_mag = star_g_int - 5*np.log10(distance_pc/10) #distance  modulus calculation for absolute mag

            
            #determined a suitable cut off for these values to get good measurements
            if ruwe < 1.1 and phot_excess <((0.15*star_colour)+1.2):
                
                name.append(data[0])
                colour.append(star_colour) #add to list of all colours
                distance.append(distance_pc)
                abs_mag_g.append(abs_mag) #add to list
                data = [float(i) for i in data]
                table.append(data)
                
            else:
                continue         

print("stars in gaia file ", stars)
table = np.array(table) #make the lists of data for each star into an array

            
def full_HR():
    #function plots the HR diagram of all the good measurements of stars
    #simple scatter plot

    print("number of stars in graph ", len(name))
    plt.plot(colour, abs_mag_g,'ro', markersize=1)
    plt.ylim(17.5, min(abs_mag_g)) #can change x and y limits to zoom in on specific areas
    plt.xlim(-1,6)
    #plt.xlim(max(colour), min(colour))
    plt.xlabel('Colour BP-RP')
    # Set the y axis label of the current axis.
    plt.ylabel('Absolute G Magnitude')
    # Set a title of the current axes.
    plt.title('Hertzsprung Russell Diagram (ruwe = 1.1, y-intercept = 1.2)')
    
    # Display a figure.
    plt.show()



def HR(colour_list, mag_list, type):
    #same as above but this time plots a line across that indicates where the white dwarfs
    #are separated from the rest of the HR diagram
    
    label = str(type)
    
    x_coord = [-1,0,1,2,3,4,5,6]
    mag_line = [6.923,10,13.077,16.154,19.231,22.308,25.385,28.462]

    print("number of stars in {} HR diagram ".format(label), len(colour_list))
    plt.plot(colour_list, mag_list,'ro', markersize=1)
    plt.plot(x_coord, mag_line, markersize=100)
    plt.ylim(17.5,-2) #again, can adjust limits to zoom in
    plt.xlim(-1,6)
    plt.xlabel('Colour BP-RP')
    # Set the y axis label of the current axis.
    plt.ylabel('Absolute G Magnitude')
    # Set a title of the current axes.
    plt.title('{} Hertzsprung Russell Diagram'.format(label))
    # Display a figure.
    plt.show()    


def WD_filter(colour, abs_mag_b, name):
    
    #goes through the full lists of colour, magnitudes and names and selects only those
    #deemed white dwarfs, appends them to the lists defined earlier
    
    for i in range(len(colour)):

        if abs_mag_b[i] > ((3.0769*colour[i])+10):
            wd_only_mag.append(abs_mag_b[i])
            wd_only_colour.append(colour[i])
            wd_only_names.append(name[i])
            
        else:
            continue
    
    return wd_only_colour, wd_only_mag, wd_only_names



def HR_histogram(colour_list, mag_list, type):
    
    #plots the same HR diagrams as previous functions but now does it in histogram form
    #this ensures that you can see the density of stars in crowded regions
    
    label = str(type)
    
    # Creating dataset
    x = colour_list
    y = mag_list
    
    # Creating bins
    x_min = np.min(x)
    x_max = np.max(x)

  
    y_min = np.min(y)
    y_max = np.max(y)

  
    x_bins = np.linspace(x_min, x_max, 1000) #creates 1000 bins between the min and max
    y_bins = np.linspace(y_min, y_max, 1000) #can adjust 
    
    fig, ax = plt.subplots(figsize =(10, 7))
    # Creating plot
    plt.hist2d(x, y, bins =[x_bins, y_bins],norm=mcolors.PowerNorm(0.3),cmap='inferno')
    plt.title('{} Hertzsprung Russell Diagram'.format(label))
    x_coord = [-1,0,1,2,3,4,5,6] #plots the white dwarf cut off line onto the histogram
    mag_line = [6.923,10,13.077,16.154,19.231,22.308,25.385,28.462]
    plt.plot(x_coord, mag_line)
    plt.ylim(np.max(y), np.min(y))
    plt.xlim(np.min(x),np.max(x))
    plt.colorbar()

    ax.set_xlabel('Colour BP-RP') 
    ax.set_ylabel('Absolute G Magnitude') 
  
    # show plot
    plt.tight_layout() 
    plt.show()
    

def close_HR_histogram(colour_list, mag_list, type):
    #same as previous function but zoomed in to a specific area
    label = str(type)
    
    # Creating dataset
    x = colour_list
    y = mag_list
    
    # Creating bins
    x_min = np.min(x)
    x_max = np.max(x)

    y_min = np.min(y)
    y_max = np.max(y)

  
    x_bins = np.linspace(x_min, x_max, 1000)
    y_bins = np.linspace(y_min, y_max, 1000)
    
    fig, ax = plt.subplots(figsize =(10, 7))
    # Creating plot
    plt.hist2d(x, y, bins =[x_bins, y_bins],norm=mcolors.PowerNorm(0.2),cmap='inferno')
    plt.title('White Dwarf Hertzsprung Russell Diagram')
    plt.ylim(np.max(y), np.min(y))
    plt.ylim(16,10.2) #could change the zoom here again
    plt.xlim(-0.5,1.5)
    plt.colorbar()

    ax.set_xlabel('Colour BP-RP') 
    ax.set_ylabel('Absolute G Magnitude') 
  
    # show plot
    plt.tight_layout() 
    plt.show()
    
def write_data(colour_list, mag_list, names_list, type):
    #writes the data to a txt file to be processed in another python script
    label = str(type)
    
    colourmag_file = open("{}_data_g.txt".format(label), "w")
    
    for i in range(len(colour_list)):
        colourmag_file.write("{} {} {}\n".format(names_list[i], colour_list[i], mag_list[i]))
    

#creating the full HR diagram
HR(colour, abs_mag_g, 'Full')

#full histogram
HR_histogram(colour, abs_mag_g, 'Full')


print("length star names list ", len(name))



#should now filter out any non white dwarves from the list, create new list
wd_only_colour, wd_only_mag, wd_only_names = WD_filter(colour, abs_mag_g,name)


#creating the WD HR diagram
HR(wd_only_colour, wd_only_mag, 'White Dwarf')

#white dwarf histogram
HR_histogram(wd_only_colour, wd_only_mag, 'White Dwarf')

close_HR_histogram(wd_only_colour, wd_only_mag, 'White Dwarf')

#write wd data to new file so that's separated
write_data(wd_only_colour, wd_only_mag, wd_only_names, 'White Dwarf')

