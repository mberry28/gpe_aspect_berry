#!/usr/bin/env python
#Created by Michael Berry 

# import modules
import numpy as np;  import matplotlib.pyplot as plt; from os import system as sys

#Numerical model
model_depth     = 400e3
model_width     = 1200e3
plateau_height  = 5e3
sample_rate_x   = 1200
sample_rate_y   = 400

dz           = np.array([20.e3,10.e3,70.e3,300.e3])
rho          = np.array([2700,2900,3230,3700]) 
plastic      = np.array([300e3,900e3,100e3])  #left, right, elevation above

X1           = 350e3 #ramp start
X2           = 400e3 #plateau start
X3           = 800e3 #plateau stop
X4           = 850e3 #ramp end 

grid_res_y   = model_depth/sample_rate_y
grid_res_x   = model_width/sample_rate_x 
H            = plateau_height / grid_res_y   #additional grid points above reference elevation
root         = ((plateau_height*rho[0])/(rho[2]-rho[1])) /grid_res_y
am_h         = int(dz[3]/grid_res_y)

slope_plat   = plateau_height/((X2-X1))
slope_H      = H/grid_res_x
slope_root   = (root/((X2-X1)/grid_res_x))

x            = np.linspace(0,model_width,num=sample_rate_x)
y            = np.linspace(0,model_depth + plateau_height,num=sample_rate_y+H)

total_column = (model_depth + plateau_height) / grid_res_y

array1  = np.zeros((int(total_column),int(sample_rate_x)))
array2  = np.zeros((int(total_column),int(sample_rate_x)))
array3  = np.zeros((int(total_column),int(sample_rate_x)))
array4  = np.zeros((int(total_column),int(sample_rate_x)))
array5  = np.zeros((int(total_column),int(sample_rate_x))) 
arrayR  = np.zeros((int(total_column),int(sample_rate_x)))

iso_eq  = 1 #1 if want to make a crustal root in isostatic equilibrium 

for j in range(int(sample_rate_x)):
    if x[j] <= X1 or x[j] >= X4:
        uc_h  = int(round(dz[0]/grid_res_y)) 
        lc_h  = int(round(dz[1]/grid_res_y))
        ml_h  = int(round(dz[2]/grid_res_y))
        
    elif x[j] > X1 and x[j] < X2:
        uc_h  = int(round(dz[0]/grid_res_y  + (j-(X1/ grid_res_x))*slope_plat))
        if (iso_eq == 1 and uc_h > int(round(dz[0]/grid_res_y))): 
            lc_h = int(round(dz[1]/grid_res_y + (j-(X1/ grid_res_x))*slope_root))
            ml_h = int(round((dz[2]/grid_res_y) - (j-(X1/ grid_res_x))*slope_root))
        else:
            lc_h = int(round(dz[1]/grid_res_y))
            ml_h = int(round(dz[2]/grid_res_y))
        h   = uc_h+lc_h+ml_h
    
    elif x[j] >= X2 and x[j] <= X3:
        uc_h = int(round(dz[0]/grid_res_y  + H))
        if (iso_eq == 1 and uc_h > int(round(dz[0]/grid_res_y))): 
            lc_h = int(round(dz[1]/grid_res_y + root))
            ml_h = int(round((dz[2]/grid_res_y) -root))
        else:
            lc_h = int(round(dz[1]/grid_res_y))
            ml_h = int(round(dz[2]/grid_res_y)) 

    elif x[j] >  X3 and x[j]< X4:
        uc_h  = int(round(((dz[0]/grid_res_y)+H) - ((j-(X3/ grid_res_x))*slope_plat)))
        if (iso_eq == 1 and uc_h > int(round(dz[0]/grid_res_y))): 
            lc_h = int(round((dz[1]/grid_res_y+root) - (j-(X3/ grid_res_x))*slope_root))
            ml_h = int(round(((dz[2]/grid_res_y)-root) + (j-(X3/ grid_res_x))*slope_root))
    else:
        lc_h = int(round(dz[1]/grid_res_y))
        ml_h = int(round(dz[2]/grid_res_y)) 
    
    h   = uc_h+lc_h+ml_h+am_h
    
    array2[:h,j] = np.column_stack([np.zeros((1,am_h)),np.zeros((1,ml_h)),np.zeros((1,lc_h)),np.ones((1,uc_h))])
    array3[:h,j] = np.column_stack([np.zeros((1,am_h)),np.zeros((1,ml_h)),np.ones((1,lc_h)),np.zeros((1,uc_h))])
    array4[:h,j] = np.column_stack([np.zeros((1,am_h)),np.ones((1,ml_h)),np.zeros((1,lc_h)),np.zeros((1,uc_h))])
    array5[:h,j] = np.column_stack([np.ones((1,am_h)),np.zeros((1,ml_h)),np.zeros((1,lc_h)),np.zeros((1,uc_h))])
    arrayR[:h,j] = np.column_stack([np.ones((1,am_h))*rho[3],np.ones((1,ml_h))*rho[2],np.ones((1,lc_h))*rho[1],np.ones((1,uc_h))*rho[0]])
    
    for i in range(int(total_column)):
        if (x[j] > plastic[0] and x[j] < plastic[1] and y[i] < plastic[2]):
            array1[i,j] = 1
        else: 
            array1[i,j] = 0

file = open('density.dat','w') 
file.write('# Only next line is parsed in format: [nx] [ny] because of keyword "POINTS:"\n') 
file.write("# POINTS: 1200 405\n") 
file.write("# Columns: x y plastic upper_crust lower_crust lithos_mantle astheno_mantle density [kg/m3]\n")
for i in range(0,405):
           for j in range(0,1200):
               file.write("%.2f %.2f %1d %1d %1d %1d %1d %.2f\n" % 
                          (x[j],y[i],array1[i,j],array2[i,j],array3[i,j],array4[i,j],array5[i,j],arrayR[i,j]))
file.close()
