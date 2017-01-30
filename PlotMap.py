# PARA HACER MAPA DE TODO EN ECUATORIALES

#- See more at: http://pelson.github.io/from_levels_and_colors.html#sthash.1whJz5IU.dpuf
#- See more at: http://pelson.github.io/from_levels_and_colors.html#sthash.1whJz5IU.dpuf
#elapsed_time()
#Usage python PlotMap.py TA
import numpy as np
import matplotlib.pyplot as plt
from my_allskymap import AllSkyMap
import astroconvert
from astropysics.coords import ICRSCoordinates,GalacticCoordinates,SupergalacticCoordinates
import sys

################### OPTIONS ###################
Type=sys.argv[1]

granularity=10
#Type='Sim_TA'
#Type='Sim_Auger'
#Type = 'S_LM_TA'
#Type = 'S_LM_Auger'
#Type = 'P_value_TA'
#Type = 'P_value_Auger'
#Type = 'TA'
#Type = 'Auger'
#Type = 'Desviacion_GAL'
#Type = 'Desviacion_EQ'
################### OPTIONS ###################

N_lon=360
N_lat=180

print "Type=",Type,"granularity=",granularity,"N_lon=",N_lon,"N_lat=",N_lat
if (Type == 'Sim_TA'  ):
    eta=72./100000.
    data2=eta*np.loadtxt("N_off.dat")
elif (Type == 'Sim_Auger' ):
    eta=201./100000.
    data2=eta*np.loadtxt("N_off.dat")
elif (Type == 'TA' or Type == 'Auger'):  
    data2=np.loadtxt("N_on.dat")
elif (Type == 'S_LM_TA' or Type == 'S_LM_Auger'):  
    data2=np.loadtxt("S_LM.dat")
elif (Type == 'P_value_TA' or Type == 'P_value_Auger'):  
    data2=np.loadtxt("P_VALUE.dat")
elif (Type == 'Desviacion_GAL'):  
    #data2=np.loadtxt("Desviacion_GAL_p60g5.dat")
    data2=np.loadtxt("Desviacion_GAL_O60g5.dat")
    #data2=np.loadtxt("Desviacion_GAL_Fe60g5.dat")
elif (Type == 'Desviacion_EQ'):  
    #data2=np.loadtxt("Desviacion_EQ_p60g5.dat")
    #data2=np.loadtxt("Desviacion_EQ_O60g5.dat")
    data2=np.loadtxt("Desviacion_EQ_Fe60g5.dat")
    
fig = plt.figure(figsize=(12,7))
main_ax = plt.axes([0.05, .19, .9, .75])  # rect=L,B,W,H

# Set up the projection and draw a grid.
if (Type == 'Desviacion_GAL'): 
    m = AllSkyMap(ax=main_ax, projection='hammer',lat_0 = 0, lon_0 = 0)
else:
    m = AllSkyMap(ax=main_ax, projection='hammer',lat_0 = 0, lon_0 = 180)


m.drawmapboundary(fill_color='white')

m.drawparallels(np.arange(-60,61,30), linewidth=1, dashes=[1,2],
    labels=[1,0,0,0], fontsize=20,color='w')
m.drawmeridians(np.arange(0,361,60), linewidth=1, dashes=[1,2],color='w')
    
if (Type != 'Desviacion_GAL' ): 
    lons = np.arange(0,361,60)
    m.label_meridians(lons, fontsize=18, vnudge=0,halign='left', hnudge=0,color='w')  # nudge<0 shifts to right
else:
    lons1 = np.arange(-180,181,60)
    m.label_meridians(lons1, fontsize=18, vnudge=0,halign='left', hnudge=0,color='k')  # nudge<0 shifts to right
x=np.linspace(0,360,num=N_lon*granularity)
y=np.linspace(90,-90,num=N_lat*granularity)

XX, YY = m(*np.meshgrid(x-180,y))
if (Type == 'Desviacion_GAL'): 
    XX, YY = m(*np.meshgrid(x-180,y))
else:
    XX, YY = m(*np.meshgrid(x,y))

im=m.pcolormesh(XX,YY,data2,cmap='jet',shading='interp')
#if (Type == 'S_LM'): im.set_clim(vmin=-np.max(data2), vmax=np.max(data2))
if (Type == 'Sim_TA' ): im.set_clim(vmin=0, vmax=19)
if (Type == 'Sim_Auger'): im.set_clim(vmin=0, vmax=22)
if (Type != 'Desviacion_GAL' ): 
    cb = m.colorbar(im,"right", size="3%", pad='3%')
else:
    cb = m.colorbar(im,"right", size="3%", pad='6%')

#m.contour(XX,YY,Overlap_Datos_Gal)

#main_ax.set_title('Zenithal Acceptance between 55 and 80 deg')
main_ax.title.set_fontsize(15)
# PLANO G
phi_G=np.linspace(0, 2*np.pi, num=400)
theta_G=np.linspace(0, 0, num=400)
x1=np.zeros(theta_G.size);y1=np.zeros(theta_G.size)
for k in range(theta_G.size):
    x1[k],y1[k]=astroconvert.gal2eq(phi_G[k],theta_G[k])
    x1[k]=np.rad2deg(x1[k])
    y1[k]=np.rad2deg(y1[k])
x2=np.sort(x1, axis=0)
sort_indices = np.argsort(x1, axis=0)
y2=y1[sort_indices]
x3, y3 = m(x2,y2)
if (Type != 'Desviacion_GAL'):  
    m.plot(x3,y3,color='w',linewidth=3.0, linestyle="--")



# PLANO SG
phi_SG=np.linspace(0, 360, num=400,endpoint=False)
theta_SG=np.linspace(0, 0, num=400,endpoint=False)
x1=np.zeros(theta_SG.size);y1=np.zeros(theta_SG.size)
x2=np.zeros(theta_SG.size);y2=np.zeros(theta_SG.size)
x11=np.zeros(theta_SG.size);y11=np.zeros(theta_SG.size)

for k in range(theta_SG.size):
    a=SupergalacticCoordinates(phi_SG[k],theta_SG[k])
    x1[k]=a.convert(GalacticCoordinates).l.radians
    y1[k]=a.convert(GalacticCoordinates).b.radians
    x2[k],y2[k]=astroconvert.gal2eq(x1[k],y1[k])
    x2[k]=np.rad2deg(x2[k])
    y2[k]=np.rad2deg(y2[k])
    x11[k]=np.rad2deg(x1[k])
    y11[k]=np.rad2deg(y1[k])
x3=np.sort(x2, axis=0)
sort_indices = np.argsort(x2, axis=0)
y3=y2[sort_indices]
x4, y4 = m(x3,y3)

if (Type != 'Desviacion_GAL'):  
    m.plot(x4,y4,color='w',linewidth=3.0, linestyle="-")
else:
    indices=np.where(x11>=180)[0]
    x3=x11[indices]
    y3=y11[indices]
    x4, y4 = m(x3,y3)
    m.plot(x4,y4,color='k',linewidth=3.0, linestyle="-")
    indices=np.where(x11<180)[0]
    x3=x11[indices]
    y3=y11[indices]
    x31=np.sort(x3, axis=0)
    sort_indices = np.argsort(x3, axis=0)
    y31=y3[sort_indices]
    x4, y4 = m(x31,y31)
    m.plot(x4,y4,color='k',linewidth=3.0, linestyle="-")
    
# Linea de Lat_Eq_min
if (Type=="TA" or Type=="Sim_TA" or Type=="S_LM_TA"  or Type=='P_value_TA'):
    # Linea de Lat_Eq_min
    theta_max_deg=55 # Max zenith angle
    Lat_Obs=39.2969 # Latitud del TA
    lat_min=Lat_Obs-theta_max_deg
    #lat_min=-10
    x1=np.linspace(0, 360, num=400)
    y1=np.linspace(lat_min, lat_min, num=400)
    x2, y2 = m(x1,y1)
    m.plot(x2,y2,color='w',linewidth=2.0, linestyle="-.")
if (Type=="Auger" or  Type=="Sim_Auger" or  Type=="S_LM_Auger" or Type=='P_value_Auger'):
    # Linea de Lat_Eq_max
    theta_max_deg=80 # Max zenith angle
    Lat_Obs=-35.2068
    lat_min=Lat_Obs-theta_max_deg
    lat_max=Lat_Obs+theta_max_deg
    x1=np.linspace(0, 360, num=400)
    y1=np.linspace(lat_max, lat_max, num=400)
    x2, y2 = m(x1,y1)
    m.plot(x2,y2,color='w',linewidth=2.0, linestyle="-.")


if (Type != 'Desviacion_GAL' and Type != 'Desviacion_EQ'): 
    x1=np.linspace(0, 360, num=400)
    y1=np.linspace(90-lat_min, 90-lat_min, num=400)
    x2, y2 = m(x1,y1) 
    m.plot(x2,y2,color='w',linewidth=2.0, linestyle="-.")
#Test
#xx, yy = m(-60,40)
#m.plot(xx, yy, 'bo', markersize=24)
min_value=np.min(data2)
max_value=np.max(data2)
print "Min=",min_value,"Max=",max_value

plt.savefig('Figure.png')

#plt.show()
