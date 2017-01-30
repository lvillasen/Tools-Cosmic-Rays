#!/usr/env python
# -*- coding: utf-8 -*-
# Usage python Plotevents.py TA Cat p
import os
import sys
Obs=sys.argv[1]
if len(sys.argv) == 3 and sys.argv[2] == 'Cat':
	add_catalog=1
else:
	add_catalog=0

if (len(sys.argv) == 3 and sys.argv[2] != 'Cat'):
    particle = sys.argv[2] 
elif (len(sys.argv) == 4 and sys.argv[2] == 'Cat') :
    particle = sys.argv[3]
else:
    particle = ''


R_s=3

##################### Preprocesing
if Obs=='TA':
    cmd='sh Preprocessing.sh TA'
elif Obs=='Auger':
    cmd='sh Preprocessing.sh Auger'
elif Obs=='Both':
    cmd='sh Preprocessing.sh Both'
else:
    print " "
    sys.exit("The argument passed must be 'TA' or 'Auger' or 'Both'. Try again ..")
os.system(cmd)
# Creates file Data.dat with this format:
# n θ (deg) E (EeV) α (deg) δ (deg) l(deg) b (deg)

##################### Import modules

import pandas as pd
import numpy as np
import astroconvert
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math as math
from matplotlib.colors import from_levels_and_colors 
from mpl_toolkits.basemap import Basemap, addcyclic
from matplotlib.colors import Normalize, NoNorm 
from matplotlib.colorbar import ColorbarBase
import matplotlib.ticker as ticker
from my_allskymap import AllSkyMap
from astropysics.coords import ICRSCoordinates,GalacticCoordinates,SupergalacticCoordinates
import matplotlib.cm as cm
import cosmolopy.distance as cd
cosmo = {'omega_M_0':0.286, 'omega_lambda_0':0.714, 'omega_k_0':0.0, 'h':0.696}
from matplotlib.patches import Polygon
import astroconvert



##################### Functions
def angle(phi_1,theta_1,phi_2,theta_2): # Angle in degrees, input in degrees
    x_1=np.sin(np.deg2rad(theta_1))*np.cos(np.deg2rad(phi_1))
    y_1=np.sin(np.deg2rad(theta_1))*np.sin(np.deg2rad(phi_1))
    z_1=np.cos(np.deg2rad(theta_1))
    x_2=np.sin(np.deg2rad(theta_2))*np.cos(np.deg2rad(phi_2))
    y_2=np.sin(np.deg2rad(theta_2))*np.sin(np.deg2rad(phi_2))
    z_2=np.cos(np.deg2rad(theta_2))
    return np.rad2deg(math.acos(x_1*x_2+y_1*y_2+z_1*z_2))
def my_Polygon(m,RA,DEC,radio,facecolor,linewidth,alpha):
    # LV 6/Jul/2015
    incluye_Polo=0
    if (angle(RA,90-DEC,0,0)<radio ): incluye_Polo=1 # North pole
    if (angle(RA,90-DEC,0,180)<radio ): incluye_Polo=-1 # South pole
    phi_elipse,theta_elipse=Elipse(RA,DEC,radio)
    segment = np.vstack((phi_elipse,theta_elipse)) 
    threshold = 90.                                                                                      
    isplit = np.nonzero(np.abs(np.diff(segment[0])) > threshold)[0]                                                                                         
    subsegs = np.split(segment,isplit+1,axis=+1) 
    #print np.shape(subsegs)
    #plot the subsegments 
    N_subseg=0
    for seg in subsegs:
        N_subseg+=1
   
    if (incluye_Polo==0):
        n=0
        for seg in subsegs:
            n+=1
            if (N_subseg == 3 and n==1): 
                a1=seg[0];b1=seg[1];a2=a1;b2=b1
            if (N_subseg == 3 and n==3): 
                a2= np.hstack((seg[0],a1));b2= np.hstack((seg[1],b1))
            else:
                a2=seg[0];b2=seg[1]
            a=seg.shape
            if (  len(a2) <200):
                lon_lim = -1
                if (360-a2[0]<1 or 360-a2[len(a2)-1]<1): lon_lim=360
                elif (a2[0]<1 or a2[len(a2)-1]<1) : lon_lim=0
                a=np.ones(len(a2))*lon_lim
                a2= np.hstack((a,a2))
                b=np.linspace(b2[len(b2)-1],b2[0],num=len(b2))
                b2=np.hstack((b,b2))
            x,y = m(a2,b2) 
            p = Polygon(zip(x,y),facecolor=facecolor,linewidth=linewidth,alpha=alpha)
            if (N_subseg <3 or n>1):
                plt.gca().add_patch(p)
                #m.plot(x,y,color='b',markersize=1)
    if (incluye_Polo==1 or incluye_Polo==-1):
        n=0
        for seg in subsegs:
            n+=1
            if (N_subseg == 2 and n==1): 
                a1=seg[0];b1=seg[1];a2=a1;b2=b1
            if (N_subseg == 2 and n==2): 
                a2= np.hstack((seg[0],a1));b2= np.hstack((seg[1],b1))
            else:
                a2=seg[0];b2=seg[1]
            a=seg.shape
            if (  len(a2) == 200):
                a0=np.ones(len(a2)/2)*360
                a1=np.ones(len(a2)/2)*0
                a_temp= np.hstack((a0,a1))
                if (incluye_Polo==-1): a_temp=a_temp[::-1]
                b0=np.linspace(b2[len(b2)-1],90*incluye_Polo,num=len(a2)/2)
                b1=np.linspace(90*incluye_Polo,b2[0],num=len(a2)/2)
                b_temp= np.hstack((b0,b1))
                a2= np.hstack((a2,a_temp))
                b2= np.hstack((b2,b_temp))
            x,y = m(a2,b2) 
            p = Polygon(zip(x,y),facecolor=facecolor,linewidth=linewidth,alpha=alpha)
            if (N_subseg <2 or n>1):
                plt.gca().add_patch(p)
                #m.plot(x,y,color='b',markersize=1)
    return()
def Elipse(Lon_deg,Lat_deg,radio_deg):
    # LV 6/Dic/2014
    import numpy as np
    N=200
    Theta_0_deg=Lat_deg-90.
    Theta_0=np.deg2rad(Theta_0_deg)
    Theta_eq=np.deg2rad(90-radio_deg)
    N=200
    theta_eq=Theta_eq*np.ones(N)
    theta=np.pi/2-theta_eq
    phi=np.linspace(-np.pi,np.pi,N)
    x1=np.sin(theta)*np.cos(phi)
    y1=np.sin(theta)*np.sin(phi)
    z1=np.cos(theta)
    x1_p=x1*np.cos(Theta_0)-z1*np.sin(Theta_0)
    theta_p_calc=np.arccos(np.sin(theta)*np.cos(phi)*np.sin(Theta_0) + np.cos(theta)*np.cos(Theta_0))
    phi_p_calc=np.rad2deg(np.arctan((np.sin(theta)*np.sin(phi))/(np.sin(theta)*np.cos(phi)*np.cos(Theta_0) - np.cos(theta)*np.sin(Theta_0))))
    my_phi=phi_p_calc
    msk=(x1_p<0)
    A=msk*180
    phi_p_calc=phi_p_calc+A
    msk1=(phi_p_calc>180)
    A1=msk1*360
    phi_p_calc=phi_p_calc-A1
    phi_final_deg=phi_p_calc + Lon_deg  
    for i in range(np.size(phi_final_deg)):
        phi_final_deg[i]=phi_final_deg[i]%360.
    theta_final_rad=np.pi/2-theta_p_calc
    theta_final_deg=np.rad2deg(theta_final_rad)
    return (phi_final_deg,theta_final_deg)




##################### Convert coordinates

# Converts coordinates to galactic coordinates
#n θ (deg) E (EeV) α (deg) δ (deg)
#n θ (deg) E (EeV) α (deg) δ (deg) l(deg) b (deg)

Data = np.loadtxt ("Data.dat");
Zenith_Obs=Data[:,1]  
E_Obs=Data[:,2]  
RA_Obs=Data[:,3]
DEC_Obs=Data[:,4]
GAL_l_Obs=np.zeros(RA_Obs.size)
GAL_b_Obs=np.zeros(DEC_Obs.size)

for k in range(RA_Obs.size):
    RA=RA_Obs[k]
    DEC=DEC_Obs[k]
    x2,y2=astroconvert.eq2gal(np.deg2rad(RA),np.deg2rad(DEC))
    GAL_l_Obs[k]=np.rad2deg(x2)
    GAL_b_Obs[k]=np.rad2deg(y2)
np.savetxt("GAL_l_Obs", GAL_l_Obs,fmt='%.2f')
np.savetxt("GAL_b_Obs", GAL_b_Obs,fmt='%.2f')
cmd='paste Data.dat GAL_l_Obs GAL_b_Obs > Data.dat.GAL'
os.system(cmd)
cmd='echo "n Theta(º) E(EeV) RA(º) DEC(º) l(º) b(º)"'
os.system(cmd)
cmd='head Data.dat.GAL'
os.system(cmd)
cmd='echo "Number of events";'
os.system(cmd)
cmd='wc  Data.dat | awk "{print $1}"'
os.system(cmd)

##################### Load data

if particle == '':
    Data = np.loadtxt ("Data.dat.GAL");
    backpropagated = 0
elif particle == 'p':
    Data = np.loadtxt ("Data.p.jf2012.dat");
    backpropagated = 1
elif particle == 'O':
    Data = np.loadtxt ("Data.O.jf2012.dat");
    backpropagated = 1
elif particle == 'Fe':
    Data = np.loadtxt ("Data.Fe.jf2012.dat");
    backpropagated = 1


##################### Cuts

E_max=700
if (Obs=="TA"):
    E_min=57
    theta_max_deg=55 # Max zenith angle
    Lat_Obs=39.2969 # Latitud del TA
if (Obs=="Auger" or Obs=="Both" ):
    E_min=54
    theta_max_deg=80 # Max zenith angle
    Lat_Obs=-35.2068 # Latitud del Auger
theta_min_deg=0 # Min zenith angle
if backpropagated == 0:
    l_Obs=Data[:,5]*((Data[:,2]>=E_min) & (Data[:,2]<E_max) & (Data[:,1]<theta_max_deg ) & (Data[:,1]>=theta_min_deg ))*1; 
    l_Obs=l_Obs[l_Obs!=0]
    b_Obs=Data[:,6]*((Data[:,2]>=E_min) & (Data[:,2]<E_max) & (Data[:,1]<theta_max_deg ) & (Data[:,1]>=theta_min_deg ))*1; 
    b_Obs=b_Obs[b_Obs!=0]
elif backpropagated == 1:
    l_Obs=Data[:,7]*((Data[:,2]>=E_min) & (Data[:,2]<E_max) & (Data[:,1]<theta_max_deg ) & (Data[:,1]>=theta_min_deg ))*1; 
    l_Obs=l_Obs[l_Obs!=0]
    b_Obs=Data[:,8]*((Data[:,2]>=E_min) & (Data[:,2]<E_max) & (Data[:,1]<theta_max_deg ) & (Data[:,1]>=theta_min_deg ))*1; 
    b_Obs=b_Obs[b_Obs!=0]

E_Obs=Data[:,2]*((Data[:,2]>=E_min) & (Data[:,2]<E_max) & (Data[:,1]<theta_max_deg ) & (Data[:,1]>=theta_min_deg ))*1; 
E_Obs=E_Obs[E_Obs!=0]
Zenith_Obs=Data[:,1]*((Data[:,2]>=E_min) & (Data[:,2]<E_max) & (Data[:,1]<theta_max_deg ) & (Data[:,1]>=theta_min_deg ))*1; 
Zenith_Obs=Zenith_Obs[Zenith_Obs!=0]

RA_Obs=np.zeros(l_Obs.size)
DEC_Obs=np.zeros(l_Obs.size)

for k in range(l_Obs.size):
    l=l_Obs[k]
    b=b_Obs[k]
    x2,y2=astroconvert.gal2eq(np.deg2rad(l),np.deg2rad(b))
    RA_Obs[k]=np.rad2deg(x2)
    DEC_Obs[k]=np.rad2deg(y2)
Zenith_Datos=Zenith_Obs
E_Datos=E_Obs
RA_Datos=RA_Obs
DEC_Datos=DEC_Obs
print DEC_Datos.size, "events with for E>=", E_min, "EeV and zenith <", theta_max_deg, "degree"
np.savetxt("Zenith_Datos", Zenith_Datos,fmt='%.2f')
np.savetxt("E_Datos", E_Datos,fmt='%.2f')    
np.savetxt("RA_Datos", RA_Datos,fmt='%.2f')
np.savetxt("DEC_Datos", DEC_Datos,fmt='%.2f')


##################### Plot events

if add_catalog == 1:
    RA_AGN = np.loadtxt ("RA_AGN.dat");         
    DEC_AGN = np.loadtxt ("DEC_AGN.dat"); 
    distance = np.loadtxt ("distance.dat"); 
    flux = np.loadtxt ("flux.dat"); 
    name=[]
    with open('name.dat') as f:
        for line in f:
            name.append(line.strip())            
################### OPTIONS ###################


N_lon=360
N_lat=180

fig = plt.figure(figsize=(12,7))
ax = plt.axes([0.05, .19, .9, .75])  # rect=L,B,W,H

# Set up the projection and draw a grid.
m = AllSkyMap(ax=ax, projection='hammer',lat_0 = 0, lon_0 = 180)
m.drawmapboundary(fill_color='white')
#m.drawmapboundary(fill_color='b')

m.drawparallels(np.arange(-60,61,30), linewidth=0.5, dashes=[1,2],
    labels=[1,0,0,0], fontsize=13)
m.drawmeridians(np.arange(0,361,60), linewidth=0.5, dashes=[1,2])
lons = np.arange(0,361,60)
m.label_meridians(lons, fontsize=13, vnudge=0,halign='left', hnudge=0)  # nudge<0 shifts to right
N_Datos=np.size(RA_Datos)
x1=np.zeros(N_Datos)
n=0
for k in range(np.size(RA_Datos)):
    n=n+1
    RA=RA_Datos[k]
    DEC=DEC_Datos[k]
    my_Polygon(m,RA,DEC,R_s,'r',1,.8)
#x3, y3 = m(RA_Datos,DEC_Datos)
#m.scatter(x3, y3, s=80, marker='o', linewidths=1, edgecolors='k',facecolors='b') 
if add_catalog ==1:
    R_catalog=[]
    print "Name","\t","\t","RA","\t", "DEC","\t","Distance" ,"\t","Radius"

    for k in range(RA_AGN.size):
        #R=10*16/distance[k]**2*flux[k]/100
        R=flux[k]/50
        R2=float("{0:.2f}".format(float(R)))
        R_catalog.append(R2)
        if k<15:
            print name[k],"\t","\t",RA_AGN[k],"\t", DEC_AGN[k],"\t",distance[k],"\t", R_catalog[k]
        my_Polygon(m,RA_AGN[k],DEC_AGN[k],R_catalog[k],'g',1,.8)
#ax.set_title('Equatorial Coordinates')
ax.title.set_fontsize(15)

# PLANO G
phi_G=np.linspace(0, 2*np.pi, num=200)
theta_G=np.linspace(0, 0, num=200)
x1=np.zeros(theta_G.size);y1=np.zeros(theta_G.size)
for k in range(theta_G.size):
    x1[k],y1[k]=astroconvert.gal2eq(phi_G[k],theta_G[k])
    x1[k]=np.rad2deg(x1[k])
    y1[k]=np.rad2deg(y1[k])
x2=np.sort(x1, axis=0)
sort_indices = np.argsort(x1, axis=0)
y2=y1[sort_indices]
x3, y3 = m(x2,y2)
m.plot(x3,y3,color='b',linewidth=2.0, linestyle="--")

# PLANO SG
phi_SG=np.linspace(0, 360, num=100,endpoint=False)
theta_SG=np.linspace(0, 0, num=100,endpoint=False)
x1=np.zeros(theta_SG.size);y1=np.zeros(theta_SG.size)
x2=np.zeros(theta_SG.size);y2=np.zeros(theta_SG.size)
for k in range(theta_SG.size):
    a=SupergalacticCoordinates(phi_SG[k],theta_SG[k])
    x1[k]=a.convert(GalacticCoordinates).l.radians
    y1[k]=a.convert(GalacticCoordinates).b.radians
    x2[k],y2[k]=astroconvert.gal2eq(x1[k],y1[k])
    x2[k]=np.rad2deg(x2[k])
    y2[k]=np.rad2deg(y2[k])
x3=np.sort(x2, axis=0)
sort_indices = np.argsort(x2, axis=0)
y3=y2[sort_indices]
x4, y4 = m(x3,y3)
m.plot(x4,y4,color='black',linewidth=2.0, linestyle="-")

if (Obs=="TA"):
    # Linea de Lat_Eq_min
    lat_min=Lat_Obs-theta_max_deg
    #lat_min=-10
    x1=np.linspace(0, 360, num=400)
    y1=np.linspace(lat_min, lat_min, num=400)
    x2, y2 = m(x1,y1)
    m.plot(x2,y2,color='black',linewidth=2.0, linestyle="-.")
if (Obs=="Auger"):
    # Linea de Lat_Eq_max
    lat_max=Lat_Obs+theta_max_deg
    x1=np.linspace(0, 360, num=400)
    y1=np.linspace(lat_max, lat_max, num=400)
    x2, y2 = m(x1,y1)
    m.plot(x2,y2,color='black',linewidth=2.0, linestyle="-.")

plt.savefig('Figure.png')
#plt.show()


np.savetxt("RA_Datos.dat", RA_Datos)
np.savetxt("DEC_Datos.dat", DEC_Datos)
