#!/usr/bin/env python

import numpy as np;  
import matplotlib.pyplot as plt; 
from os import system as sys
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Code for generating a ASPECT compatible geotherm with changing elevation- a plateau with ramps.
# define the points at which the ramp starts and stops.

#             _ _ _ _ _             _
#            /         \            |
#           /           \           plateau height
#_ _ _ _ _ /             \ _ _ _ _ _|
#         ^   ^       ^  ^ 
#     ramp  plat     plat ramp
#     start start    end  end
# 
# the model uses an inputted surface heatflow Q[0], surface temp T[0], thickness of the upper crust dz[0],  
# lower crust dz[1], upper mantle dz[2], surface radiogenic heat production A[0], and the temperature at which the 
# adiabatic mantle takes over T[4]. 
# 
# the thickness of the plateau occurs in the upper crust for this code.

#output file name
Filename        = 'geotherm_out.dat'

#Numerical model
model_depth     = 400e3
model_width     = 1200e3
plateau_height  = 5e3
sample_rate_x   = 1200
sample_rate_y   = 400

#Lateral ramp dimensions
X1              = 350e3 #ramp start
X2              = 400e3 #plateau start
X3              = 800e3 #plateau stop
X4              = 850e3 #ramp end 

dz              = np.array([20.e3,10.e3,5.e3,75.e3,100.e3]) #layer thicknesses
A               = np.array([1.4e-6,0,0,0])              #heat production
T               = np.array([273.,0.,0.,0.,1573.])       #T[0] is top T[3] is base
Q               = np.array([60e-3,0.,0.,0.,0.0])        #heat flow 
k               = 3.0     # Thermal conductivity W/(m*K)
ATG             = 0.3     # Asthenospheric thermal gradient (k/km) 

Q_plat          = 65e-3   #heatflow at top of plateau
Q_ref           = 60e-3   #heatflow reference 
uc_ref          = 20e3    #upper crustal reference thicknessnes

#plotting? 'Y' or 'N'
plot            = 'N'   

#------------------------------------------------End of inputs
grid_res_y      = model_depth/sample_rate_y
grid_res_x      = model_width/sample_rate_x 
H               = plateau_height / grid_res_y  # additional grid points above reference elevation
slope_plat      = plateau_height/((X2-X1))
slope_Q         = (Q_plat-Q_ref)/(X2-X1) 
x               = np.linspace(0,model_width,num=sample_rate_x)

for j in range(sample_rate_x):
    if x[j] <= X1 or x[j] >= X4:
        y     = model_depth
        dz[0] = uc_ref
        h     = sample_rate_y
        Q[0]  = Q_ref
    elif x[j] > X1 and x[j] < X2:
        dz[0] = uc_ref + ((x[j]-X1)*slope_plat)
        y     = model_depth + ((x[j]-X1)*slope_plat)
        h     = int(round(sample_rate_y + (j-(X1/grid_res_x))*slope_plat))
        Q[0]  = Q_ref + (slope_Q*(x[j]-(X1)))
    elif x[j] >=X2 and x[j] <= X3:
        dz[0] = uc_ref + (plateau_height)
        y     = model_depth + plateau_height
        h     = int(sample_rate_y + H)
        Q[0]  = Q_plat
    elif x[j] >  X3 and x[j]< X4:
        dz[0] = (uc_ref + plateau_height)- ((x[j]-X3)*slope_plat)
        y     = (model_depth + plateau_height) - ((x[j]-X3)*slope_plat)
        h     =  int(round((sample_rate_y + H) -(j-(X3/ grid_res_x))*slope_plat))
        Q[0]  = Q_plat - (slope_Q*(x[j]-(X3))) 
    
    # Calculate heat flow at the top of the upper crust
    Q[1] = Q[0] - (A[0]*dz[0])
                
    # Calculate radiogenic heat concentration at the top of the upper crust
    T[1] = T[0] + A[0]*((dz[0]**2)/(2.*k))+ ((Q[1]/k)*dz[0]) 
    
    Q[2] = Q[1] - (A[1]*dz[1])
    
    T[2] = T[1] + A[1]*((dz[1]**2)/(2.*k))+ ((Q[2]/k)*dz[1]) 
                
    Q[3] = Q[2] - (A[2]*dz[2])
                
    T[3] = T[2] + A[2]*((dz[2]**2)/(2.*k))+ ((Q[3]/k)*dz[2]) 

    
    z    = np.linspace(0,y,num=h).reshape(h,1)
    t    = np.zeros((h,1))
    X    = np.ones((int(sample_rate_y+H),1))*x[j]
    Z    = (np.linspace(0,int(sample_rate_y+H),num=int(sample_rate_y+H))*grid_res_y).reshape(int(sample_rate_y+H),1)
    tt   = np.ones((int(sample_rate_y+H),1))*273
    
    for i in range(h):
            if z[i] <=dz [0]: 
                t[i] = T[0] + (Q[0]/k)*z[i] - (A[0]*(z[i]**2))/(2*k)
            elif z[i] > dz[0] and z[i] <= (dz[0]+dz[1]):
                t[i] = T[1] + (Q[1]/k)*(z[i] - dz[0]) - (A[1]*((z[i]-dz[0])**2))/(2*k)
            elif z[i] > (dz[0]+dz[1]) and z[i] <= (dz[0]+dz[1]+dz[2]):
                t[i] = T[2] + (Q[2]/k)*(z[i]-dz[0]-dz[1]) - (A[2]*((z[i]-dz[0]-dz[1])**2))/(2*k)
            elif z[i] > (dz[0]+dz[1]+dz[2]): #and z[i]<=(dz[0]+dz[1]+dz[2]+dz[3]):
                t[i] = T[3] + (Q[3]/k)*(z[i]-dz[0]-dz[1]-dz[2]) - (A[3]*((z[i]-dz[0]-dz[1]-dz[2])**2))/(2*k)
            if t[i] >= T[4]:
                  t[i] = T[4] + ATG*(z[i]-dz[0]-dz[1]-dz[2]-dz[3])/1000

    tt[-h:]  = t  
    Z[:h]    = z
    Z        = np.flipud(Z)
    if j == 0:
        Lat  = (X)
        Ele  = Z
        Temp = tt
    else:
        Lat  = np.hstack((Lat,X))
        Ele  = np.hstack((Ele,Z))
        Temp = np.hstack((Temp,tt))
Ele  = np.flipud(Ele) 
Temp = np.flipud(Temp)

if plot == 'Y':
    plt.figure(figsize=(10,15)) 
    ax = plt.gca()
    im = ax.imshow(np.flipud(Temp),aspect=1,cmap='seismic') 
    plt.xlabel('Distance [km]')
    plt.ylabel('Depth [km]') 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05) 
    cbar = plt.colorbar(im,cax=cax) 
    cbar.set_label('Temperature [K]', rotation=270,labelpad=14)
    cbar.ax.invert_yaxis() 
    plt.savefig('geotherm_crosssection.png')
    plt.show()  

    plt.figure(figsize=(5,5))
    plt.plot(Temp[:,0],(Ele[:,0]),color='red',lw=2) 
    plt.savefig('single_geotherm.png')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Elevation [m]')
    plt.savefig('geotherm_single.jpg')
    
#plt.show()  


file = open(Filename,'w') 
file.write('# Only next line is parsed in format: [nx] [ny] because of keyword "POINTS:"\n') 
file.write("# POINTS: %4d %4d\n" %(sample_rate_x, sample_rate_y+H)) 
file.write("# Columns: x y temperature [K]\n")
for i in range(0,int(sample_rate_y+H)):
           for j in range(0,int(sample_rate_x)):
              file.write("%.2f %.2f %.2f\n" % (Lat[i,j], Ele[i,j], Temp[i,j]))
file.close()      
