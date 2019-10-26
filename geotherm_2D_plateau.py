
import numpy as np;  import matplotlib.pyplot as plt; from os import system as sys
import pandas as pd;

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


model_width     = 1200e3
sampling_points = 120 
ref_y_points    = 400
grid_resolution = 1.e3 
column_points   = 400   # grid points in the Y direction of the geotherm


# convert all values to km in advance
k      = 3.0   #Thermal conductivity W/(m*K)
ATG    = 0.3   #Asthenospheric thermal gradient (k/km)
Q_plat = 65e-3 #heatflow at top of plateau
Q_ref  = 60e-3 #heatflow reference

x      = np.linspace(0,model_width,num=sampling_points)
dz     = np.array([20.e3,10.e3,5.e3,75.e3,100.e3]) 
A      = np.array([1.4e-6,0,0,0]) #heat production
T      = np.array([273.,0.,0.,0.,1573.]) #temp
Q      = np.array([60e-3,0.,0.,0.,0.0]) #heat flow 
k      = 3.0     # Thermal conductivity W/(m*K)
ATG    = 0.3  # Asthenospheric thermal gradient (k/km)
 

H = plateau_height / grid_resolution  #plateau height -> additional grid points above reference elevation
X1 = 350 #ramp start
X2 = 400 #plateau start
X3 = 800 #plateau stop
X4 = 850 #ramp end 

slope_ramp = (X2-X1)/H
slope_Q    = (X1-X1)/Q_max-

for j in range(sampling_points)):
    if x[j] <= X1*1e3 or x[j] >= X4*1e3:
     y = 400e3
     dz[0] = 20.e3
     h = column_points 
     Q[0] = Q_ref
    elif x[j] > X1*1e3 and x[j] < X2*1e3:
     dz[0] = 20.e3+(0.1*(x[j]-(X1*1e3))) #JBN - what are terms 0.1 and 0.001? Related to ramp; I would define these variables up top and if possible calculate the values in advanced using the plateau width
     y = 400e3 + (0.01*(x[j]-(X1*1e3)))
     h = 400 + (j-35) # what is 35 related to? 35 comes from x1/number_sampling points -> calculate this in advance
     Q[0] = Q_ref + (slope_Q*(x[j]-(X1*1e3)))
    elif x[j] >=X2*1e3 and x[j] <= X3*1e3:
     dz[0] = 20.e3 + (H*1e3)
     y = 400e3 + (H*1e3)
     h = 400 + H
     Q[0] = Q_plat
    elif x[j] >  X3*1e3 and x[j]< X4*1e3:
     dz[0] = 20.e3-(slope_ramp*(x[j]-(X4*1e3)))
     y = 400e3 - (slope_ramp*(x[j]-(X4*1e3)))
     h = (400 + H) - (j-80)
     Q[0] = Q_ref - (slope_Q*(x[j]-(X4*1e3)))

    
    # Calculate heat flow at the top of the upper crust
    Q[1] = Q[0] - (A[0]*dz[0])
                
    # Calculate radiogenic heat concentration at the top of the upper crust
    T[1] = T[0] + A[0]*((dz[0]**2)/(2.*k))+ ((Q[1]/k)*dz[0]) 
    
    Q[2] = Q[1] - (A[1]*dz[1])
    
    T[2] = T[1] + A[1]*((dz[1]**2)/(2.*k))+ ((Q[2]/k)*dz[1]) 
                
    Q[3] = Q[2] - (A[2]*dz[2])
                
    T[3] = T[2] + A[2]*((dz[2]**2)/(2.*k))+ ((Q[3]/k)*dz[2]) 

    
    z = np.linspace(0,y,num=h).reshape(h,1)
    t = np.zeros((h,1))
    X = np.ones((400+H,1))*x[j]
    Z = (np.linspace(0,400+H,num=400+H)*1e3).reshape(400+H,1)
    tt = np.ones((400+H,1))*273
    
    for i in range(h):
            if z[i]<=dz[0]: 
             t[i] = T[0] + (Q[0]/k)*z[i] - (A[0]*(z[i]**2))/(2*k)
            elif z[i]>dz[0] and z[i]<=(dz[0]+dz[1]):
             t[i] = T[1] + (Q[1]/k)*(z[i]-dz[0]) - (A[1]*((z[i]-dz[0])**2))/(2*k)
            elif z[i]>(dz[0]+dz[1]) and z[i]<=(dz[0]+dz[1]+dz[2]):
             t[i] = T[2] + (Q[2]/k)*(z[i]-dz[0]-dz[1]) - (A[2]*((z[i]-dz[0]-dz[1])**2))/(2*k)
            elif z[i]>(dz[0]+dz[1]+dz[2]): #and z[i]<=(dz[0]+dz[1]+dz[2]+dz[3]):
             t[i] = T[3] + (Q[3]/k)*(z[i]-dz[0]-dz[1]-dz[2]) - (A[3]*((z[i]-dz[0]-dz[1]-dz[2])**2))/(2*k)
             if t[i] >= T[4]:
                  t[i] = T[4] + ATG*(z[i]-dz[0]-dz[1]-dz[2]-dz[3])/1000

    tt[-h:] = t  
    Z[:h]=z
    Z = np.flipud(Z)
    if j == 0:
        Lat = (X)
        Ele = Z
        Temp = tt
    else:
        Lat = np.hstack((Lat,X))
        Ele = np.hstack((Ele,Z))
        Temp = np.hstack((Temp,tt))
Ele = np.flipud(Ele) 
Temp = np.flipud(Temp)

#plt.figure(figsize=(20,20))
#plt.imshow(np.flipud(Temp),aspect=.100,cmap='seismic')

file = open('geotherm_adjusted.dat','w') 
file.write('# Only next line is parsed in format: [nx] [ny] because of keyword "POINTS:"\n') 
file.write("# POINTS: 120 405\n") 
file.write("# Columns: x y temperature [K]\n")
for i in range(0,405):
           for j in range(0,120):
              file.write("%.2f %.2f %.2f\n" % (Lat[i,j], Ele[i,j], Temp[i,j]))
file.close()      
