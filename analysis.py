#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 11:16:33 2022

@author: emilyroberts
"""

import pandas as pd
import numpy as np
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


def magindex(j):
    #converts the python index back into the actual magnitude value of a bin
    mag = (-0.25*j)+18
    
    return mag
    
def colourindex(i):
    #converts the python index back into the actual colour value of a bin
    colour = (0.125*i)-0.75
    
    return colour
    
    
#creating panda array for the data using the white dwarf txt file created from the other script
#has all the gaia data as headers for easy access
headers = ["ID", "ra", "error ra", "dec", "error dec", "parallax", "error parallax", "pm ra", "error pm ra", "pm dec", "error pm dec", "apparent g", "apparent blue", "apparent red", "ruwe", "phot excess", "ignore"]

data = np.loadtxt("whitedwarfs.txt", delimiter=" ",dtype='str')
data[data=='']=np.NaN
data2 = data
data = data.astype(np.float) #this is trying to ensure all data is same type, all columns full and equal etc
df = pd.DataFrame(data, columns=headers)
print("HERE")
print(df)


df["colour"] = df["apparent blue"]-df["apparent red"] #calculating colour

df["magnitude"] = df["apparent g"]-(5*np.log10(100/df["parallax"])) #calculating magnitude

df["ID2"] = data2[:,0].astype(str) #makes ID a string not int

print(df)



df.to_csv('pls.csv', index=False)

full_array = df.to_numpy()



items_in_row = np.size(full_array, axis = 0)
print("items in row ", items_in_row) #checking if it's got all the rows of data I need



bin_size_colour = 0.125 #these can be changed manually for finer or broader bins
bin_size_mag = 0.25     #same with magnitude bin size

colour_list = full_array[:,17] #selecting only the colour and mag columns
mag_list = full_array[:,18]

#will need to start bins here from (-0.75, 18) as that's the bottom left corner of the HR diagram
print("max colour list ", np.max(colour_list))
print("min colour list ", np.min(colour_list))
print("max mag list ", np.max(mag_list))
print("Min mag list ", np.min(mag_list))
num_colour_bins =  ((np.max(colour_list)-(-0.75))/bin_size_colour)
num_mag_bins =  ((18-np.min(mag_list))/bin_size_mag) #calculating the number of bins needed

print("num colour bins", num_colour_bins)
print("rounded ", math.ceil(num_colour_bins)) #need an integer number of bins
print("num mag bins ", num_mag_bins)
print("rounded ", math.ceil(num_mag_bins))

x_max = float(-0.75 + (math.ceil(num_colour_bins)*bin_size_colour))
print("x max ", x_max)

y_min = float(18 - (math.ceil(num_mag_bins)*bin_size_mag))
print("y min ", y_min)

x_bins = np.linspace(-0.75, x_max, math.ceil(num_colour_bins+1)) #list of the edges of the bins
y_bins = np.linspace(18, y_min, math.ceil(num_mag_bins)+1)


print("x_bins ", x_bins)
print("y_bins ", y_bins)

#assigns each star to a mag bin based on the value of the magnitude for that star
mbins = np.digitize(full_array[:,18], y_bins) 
mbins[mbins>num_mag_bins] = num_mag_bins
print("mbins ", mbins)

#assigns each star to a colour bin based on the value of the colour for that star
cbins = np.digitize(full_array[:,17], x_bins)
cbins[cbins>num_colour_bins] = num_colour_bins
print("cbins ", cbins)
print("max cbin ", np.max(cbins))

bin_indices = np.vstack((cbins-1,mbins-1)).T #adds the bin indices as a new column on the array
full_array = np.hstack((full_array, bin_indices))


items_in_row = np.size(full_array, axis = 0)
print("items in row ", items_in_row) #checks that array  manipulatin has worked, right no. of rows/columns


mag_indexing = [*range(0, np.max(mbins), 1)]
print("mag indexing ", mag_indexing)
colour_indexing = [*range(0, np.max(cbins), 1)]
print("colour indexing ", colour_indexing)
stars_done = 0

avg_A_list = []
avg_P_list = []

#lists for plotting
avg_V_u = []
avg_V_v = []
avg_V_w = []
mag_list = []
colour_list = []
#list for computing velocity dispersion S^2
square_p_list = []
powerfour_p_dash_list = []
#creating empty arrays for all values necessary, right sizes
s_squared_array = np.zeros((len(colour_indexing),len(mag_indexing)))
s_array = np.zeros((len(colour_indexing), len(mag_indexing)))
variance_u_array = np.zeros((len(colour_indexing), len(mag_indexing)))
variance_v_array = np.zeros((len(colour_indexing), len(mag_indexing)))
variance_w_array = np.zeros((len(colour_indexing), len(mag_indexing)))
variance_s_squared_array = np.zeros((len(colour_indexing), len(mag_indexing)))
sigma_s_array = np.zeros((len(colour_indexing), len(mag_indexing)))
avg_v_array = np.zeros((len(colour_indexing), len(mag_indexing)))
avg_u_array = np.zeros((len(colour_indexing), len(mag_indexing)))
avg_w_array = np.zeros((len(colour_indexing), len(mag_indexing)))

no_stars_list = 0
not_enough_stars =0 #numbers for tracking progress
few_stars = 0
bin_no = 0


check_number = []
check_vd = []

position_l_q = [] #lists for collecting the positions of the anomolous points of data
position_b_q = []

for i in range(len(colour_indexing)):
    bin_colour = full_array[full_array[:,20]==i]
    print("bin colour ", bin_colour[:,20])
    for j in range(len(mag_indexing)):
        bin_colour_mag = bin_colour[bin_colour[:,21]==j] #filters to get one magnitude-colour combo

        print(" ")
        print("  ")
        print("  ")
        print("NEW BIN")
        print("i've selected a bin with mag ", magindex(j), " and colour ", colourindex(i))
        bin_no += 1
        print("this is bin number ", bin_no)

            
        
        print("len bin colour mag ", len(bin_colour_mag))

        if len(bin_colour_mag)==0:
            print("no stars in this bin")
            no_stars_list +=1
        elif len(bin_colour_mag)==1:
            print("not enough stars in this bin")
            not_enough_stars += 1
            
        if len(bin_colour_mag)<10:
            print("fewer than ten stars in this bin")
            few_stars += 1
        else:
            A_list = []
            P_list = []
            ID_array = bin_colour_mag[:,0] #finding the column for each piece of data
            ra_array = bin_colour_mag[:,1]
            dec_array = bin_colour_mag[:,3]
            parallax_array = bin_colour_mag[:,5]
            pm_ra = bin_colour_mag[:,7]
            pm_dec = bin_colour_mag[:,9]
            print("number of stars in this bin ", len(ID_array))
            stars_done += int(len(ID_array)) #keeping track of progress
            print(stars_done)
            for k in range(len(ID_array)):
                
                parallax_as = parallax_array[k]/1000
                distance_pc = 1/parallax_as
                c = SkyCoord(ra=ra_array[k]*u.degree, dec = dec_array[k]*u.degree, pm_ra_cosdec=pm_ra[k]*u.mas/u.yr, pm_dec=pm_dec[k]*u.mas/u.yr, radial_velocity=0*u.km/u.s, distance=distance_pc*u.pc, frame='icrs')
                gc = c.galactic #converting position coords to galactic coords
                if magindex(j)==12.75: #investigating the anomolous data points in these bins
                    if colourindex(i)==-0.5:
                        position_l_q.append(gc.l.value)
                        position_b_q.append(gc.b.value)
                    if colourindex(i)==-0.25:
                        position_l_q.append(gc.l.value)
                        position_b_q.append(gc.b.value)
                        
                #this is implementing eqn 2 from the dehnen and binney 1998 paper              
                pm_l_rpy = (gc.pm_l_cosb.value)*(math.pi)*(1/648000000) #converting from mas/yr to rads/yr
                pm_b_rpy = (gc.pm_b.value)*(math.pi)*(1/648000000) #converting from mas/yr to rads/yr

                
                pmv1 = -(np.sin(math.radians(gc.l.value)))*(pm_l_rpy) - (np.cos(math.radians(gc.l.value)))*(np.sin(math.radians(gc.b.value)))*(pm_b_rpy)
                pmv2 = (np.cos(math.radians(gc.l.value)))*(pm_l_rpy) - (np.sin(math.radians(gc.l.value)))*(np.sin(math.radians(gc.b.value)))*(pm_b_rpy)
                pmv3 = (np.cos(math.radians(gc.b.value)))*(pm_b_rpy)
                #convert to radians per year   need p in km per s ultimately                           
                proper_motion_velocity = (1/parallax_as)*np.array((pmv1, pmv2, pmv3))*(977799.3256774914) #convert from pc/yr to km/s
                #print("pmv for the star ")
                #print(proper_motion_velocity)
                P_list.append(proper_motion_velocity)
                
                #finding the unit vector to star
                position = gc.cartesian
                #print("distance pc ", distance_pc)
                #print("position ", position)
                r_vector = np.array([position.x.value, position.y.value,position.z.value])
                r_unit_vector = (1/distance_pc)*r_vector
                #print("unit vector ", r_unit_vector) #in parsecs
                
                #calculating A for each star
                #equation 4 in the paper
                A = np.identity(3) - np.outer(r_unit_vector, r_unit_vector)
                #print("A ", A)
                A_list.append(A)

                
            #print("For this bin, the average A is ")
            avg_A = np.mean(np.stack(A_list), axis=0)
            #print(avg_A)
            avg_A_list.append(avg_A)
            
            #print("and the average p is ")
            avg_P = np.mean(np.stack(P_list), axis=0)
            #print(avg_P)
            avg_P_list.append(avg_P)

            #equation 5 from paper, inverting the A matrix to get average v
            invert_avg_A = np.linalg.inv(avg_A)
            
            avg_v = np.matmul(invert_avg_A, avg_P)
            print("avg_v is ", avg_v, "for the bin of colour ", i, "and mag ", j)
            #need these for plotting V components against colour/mag
            avg_V_u.append(avg_v[0])
            avg_u_array[i,j] = avg_v[0]
            avg_V_v.append(avg_v[1])
            avg_v_array[i,j] = avg_v[1]
            avg_V_w.append(avg_v[2])
            avg_w_array[i,j] = avg_v[2]
            mag_list.append(magindex(j))
            colour_list.append(colourindex(i))
            
            #equation 6 from paper, finding motion relative to the mean
            #not doing it for v because we don't know the true v, only doing it for p            
            for k in range(len(ID_array)):
                relative_motion_p = P_list[k] - np.matmul(A_list[k], avg_v)
                #print("relative motion for star ID ", int(ID_array[k]))
                #print(relative_motion_p)
                square_p_dash = np.dot(relative_motion_p, relative_motion_p)
                powerfour_p_dash = square_p_dash**2
                square_p_list.append(square_p_dash)
                powerfour_p_dash_list.append(powerfour_p_dash)
            
            s_squared = np.mean(np.stack(square_p_list),axis=0)
            powerfour_part = np.mean(np.stack(powerfour_p_dash_list), axis=0)
            print("S squared for bin colour ", i, "and mag ", j, " is ", s_squared)
            print("power four part is ", powerfour_part)
            print("number of stars is ",len(ID_array))
            s_squared_array[i,j] = s_squared
            s = np.sqrt(s_squared)
            s_array[i,j]= s   
            if magindex(j)==12.75:#also looking at anomolous data points
                # if -0.5<colourindex(i)<-0.376:
                check_number.append(len(ID_array))
                check_vd.append(s)
            
            #variance of <v>
            variance_avg_v = (1/len(ID_array))*invert_avg_A*s_squared
            #print("variance of <v> ", variance_avg_v)
            #variance_v_array[i,j] = variance_avg_v
            variance_u_array[i,j] = np.sqrt(variance_avg_v[0,0])
            variance_v_array[i,j] = np.sqrt(variance_avg_v[1,1])
            variance_w_array[i,j] = np.sqrt(variance_avg_v[2,2])
            #variance of s squared
            variance_s_squared = (1/len(ID_array))*(powerfour_part - (s_squared**2))
            print("variance of S^2 ", variance_s_squared)
            variance_s_squared_array[i,j] = variance_s_squared
            sigma_s = (1/(2*s))*np.sqrt(variance_s_squared)
            sigma_s_array[i,j] = sigma_s
print("number of average A matrices ", len(avg_A_list))
print("number of avg P vectors ", len(avg_P_list))
print("number of bins with no stars in ", no_stars_list)
print("number of bins with not enough stars ", not_enough_stars)
print("number of bins with less than 10 stars ", few_stars)
print("numbe ro fstars processed ", stars_done)
print("Min value of S squared ", np.min(s_squared_array))
print("Max value of S squared ", np.max(s_squared_array))

print("min value of variance of s squared ", np.min(variance_s_squared_array))
print("max value of variance of s squared ", np.max(variance_s_squared_array))




colour_values = ['-0.5', '0','0.5','1.0','1.5','2'] #list of strings to set the axes ticks
mag_values = ['18','16','14','12','10']


# S ARRAY
plot_array = np.transpose(s_array)
fig, ax = plt.subplots()
im = ax.imshow(plot_array, norm=mcolors.PowerNorm(1), cmap='inferno')
ax.invert_yaxis()
ax.set_xticks(np.arange(min(colour_indexing), max(colour_indexing)+1,step=4))
ax.set_yticks(np.arange(min(mag_indexing), max(mag_indexing)+1,step=8))
ax.set_xticklabels(colour_values)
ax.set_yticklabels(mag_values)
ax.set_title("Velocity Dispersion S (kms$^{-1}$)")
ax.set_xlabel("Colour BP-RP")
ax.set_ylabel("Absolute G Magnitude")
fig.colorbar(im)



#SIGMA S ARRAY
sigma_s_array = np.transpose(sigma_s_array)
fig, ax = plt.subplots()
sigma = ax.imshow(sigma_s_array,norm=mcolors.PowerNorm(0.5), cmap='inferno')
ax.invert_yaxis()
ax.set_xticks(np.arange(min(colour_indexing), max(colour_indexing)+1,step=4))
ax.set_yticks(np.arange(min(mag_indexing), max(mag_indexing)+1,step=8))
ax.set_xticklabels(colour_values)
ax.set_yticklabels(mag_values)
ax.set_title("Standard Deviation of S (kms$^{-1}$)")
ax.set_xlabel("Colour BP-RP")
ax.set_ylabel("Absolute G Magnitude")
fig.colorbar(sigma)

#V COMPONENT ARRAY
v_array = np.transpose(avg_v_array)
fig, ax = plt.subplots()
im2 = ax.imshow(v_array, norm=mcolors.PowerNorm(0.4), cmap='inferno')
ax.invert_yaxis()
ax.set_xticks(np.arange(min(colour_indexing), max(colour_indexing)+1,step=4))
ax.set_yticks(np.arange(min(mag_indexing), max(mag_indexing)+1,step=8))
ax.set_xticklabels(colour_values)
ax.set_yticklabels(mag_values)
ax.set_title("V component of <V> (kms$^{-1}$)")
ax.set_xlabel("Colour BP-RP")
ax.set_ylabel("Absolute G Magnitude")
fig.colorbar(im2)

#V STD DEV ARRAY
v_variance_array = np.transpose(variance_v_array)
fig, ax = plt.subplots()
vvar = ax.imshow(v_variance_array, norm=mcolors.PowerNorm(0.5), cmap='inferno')
ax.invert_yaxis()
ax.set_xticks(np.arange(min(colour_indexing), max(colour_indexing)+1,step=4))
ax.set_yticks(np.arange(min(mag_indexing), max(mag_indexing)+1,step=8))
ax.set_xticklabels(colour_values)
ax.set_yticklabels(mag_values)
ax.set_title("Standard Deviation of V component of <V> (kms$^{-1}$)")
ax.set_xlabel("Colour BP-RP")
ax.set_ylabel("Absolute G Magnitude")
fig.colorbar(vvar)
   
#U COMPONENT ARRAY 
u_array = np.transpose(avg_u_array)
fig, ax = plt.subplots()
im3 = ax.imshow(u_array, norm=mcolors.PowerNorm(0.5, vmin=-50, vmax=50, clip=True), cmap='inferno')
ax.invert_yaxis()
ax.set_xticks(np.arange(min(colour_indexing), max(colour_indexing)+1,step=4))
ax.set_yticks(np.arange(min(mag_indexing), max(mag_indexing)+1,step=8))
ax.set_xticklabels(colour_values)
ax.set_yticklabels(mag_values)
ax.set_title("U component of <V> (kms$^{-1}$)")
ax.set_xlabel("Colour BP-RP")
ax.set_ylabel("Absolute G Magnitude")
fig.colorbar(im3)

#U ST DEV ARRAY
u_variance_array = np.transpose(variance_u_array)
fig, ax = plt.subplots()
uvar = ax.imshow(u_variance_array, norm=mcolors.PowerNorm(0.5), cmap='inferno')
ax.invert_yaxis()
ax.set_xticks(np.arange(min(colour_indexing), max(colour_indexing)+1,step=4))
ax.set_yticks(np.arange(min(mag_indexing), max(mag_indexing)+1,step=8))
ax.set_xticklabels(colour_values)
ax.set_yticklabels(mag_values)
ax.set_title("Standard Deviation of U component of <V> (kms$^{-1}$)")
ax.set_xlabel("Colour BP-RP")
ax.set_ylabel("Absolute G Magnitude")
fig.colorbar(uvar)

#W COMPONENT ARRAY
w_array = np.transpose(avg_w_array)
fig, ax = plt.subplots()
im4 = ax.imshow(w_array, norm=mcolors.PowerNorm(0.52, vmin=-50, vmax=50), cmap='inferno')
ax.invert_yaxis()
ax.set_xticks(np.arange(min(colour_indexing), max(colour_indexing)+1,step=4))
ax.set_yticks(np.arange(min(mag_indexing), max(mag_indexing)+1,step=8))
ax.set_xticklabels(colour_values)
ax.set_yticklabels(mag_values)
ax.set_title("W component of <V> (kms$^{-1}$)")
ax.set_xlabel("Colour BP-RP")
ax.set_ylabel("Absolute G Magnitude")
fig.colorbar(im4)

#W ST DEV ARRAY
w_variance_array = np.transpose(variance_w_array)
fig, ax = plt.subplots()
wvar = ax.imshow(w_variance_array, norm=mcolors.PowerNorm(0.5), cmap='inferno')
ax.invert_yaxis()
ax.set_xticks(np.arange(min(colour_indexing), max(colour_indexing)+1,step=4))
ax.set_yticks(np.arange(min(mag_indexing), max(mag_indexing)+1,step=8))
ax.set_xticklabels(colour_values)
ax.set_yticklabels(mag_values)
ax.set_title("Standard Deviation of W component of <V> (kms$^{-1}$)")
ax.set_xlabel("Colour BP-RP")
ax.set_ylabel("Absolute G Magnitude")
fig.colorbar(wvar)


print("number stars in each bin")
print(check_number)
print("velocity dispersion")
print(check_vd)


plt.plot(position_l_q, position_b_q, 'ro', markersize=3)
plt.xlim(0,360)
plt.ylim(-90,90)
plt.title("Position of Stars in Outlying High S Bins")
plt.xlabel("$\ell$ coordinate ($\degree$)")
plt.ylabel("$b$ coordinate ($\degree$)")
plt.show()