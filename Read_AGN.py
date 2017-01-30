import pandas as pd
import numpy as np
import cosmolopy.distance as cd
import astroconvert
cosmo = {'omega_M_0':0.286, 'omega_lambda_0':0.714, 'omega_k_0':0.0, 'h':0.696}

VCV=pd.read_csv('VCV_z03.dat', sep='|', skiprows=4)
distance=[];RA_AGN=[];DEC_AGN=[];name=[];flux=[]
j=0
for i in range(len(VCV)):
    if VCV.values[i][2] >0 and VCV.values[i][2] <.03:
        RA,DEC =np.rad2deg(astroconvert.gal2eq(np.deg2rad(VCV.values[i][4]),np.deg2rad(VCV.values[i][5])) )      
        RA2=float("{0:.2f}".format(float(RA)))
        DEC2=float("{0:.2f}".format(float(DEC)))
        RA_AGN.append(RA2)
        DEC_AGN.append(DEC2)
        D=cd.comoving_distance(VCV.values[i][2], **cosmo)
        D2=float("{0:.2f}".format(float(D)))
        distance.append(D2)
        name.append(VCV.values[i][1])
        flux.append(100.0)
        j+=1
sorted_index, sorted_d = zip(*sorted([(i,e) for i,e in enumerate(distance)], key=lambda x:x[1]))
RA_AGN_sorted=[RA_AGN[i] for i in sorted_index]
DEC_AGN_sorted=[DEC_AGN[i] for i in sorted_index]
distance_sorted=[distance[i] for i in sorted_index]  
name_sorted=[name[i] for i in sorted_index]    

print  "Total number of objects=",len(VCV)
print  "Total number of objects with measured redshift < 0.03=",j

print "Name","\t","\t","\t","RA","\t","DEC","\t","\t","distance (Mpc)"
for i in range(10):
    print name[i],"\t", RA_AGN_sorted[i],"\t", DEC_AGN_sorted[i],"\t", distance_sorted[i],"\t"
print " ......."
print name[-1],"\t", RA_AGN_sorted[-1],"\t", DEC_AGN_sorted[-1],"\t",distance_sorted[-1]
np.savetxt("RA_AGN.dat", RA_AGN_sorted)
np.savetxt("DEC_AGN.dat", DEC_AGN_sorted)
np.savetxt("distance.dat", distance_sorted)
np.savetxt("flux.dat", flux)
with open ("name.dat","w")as fp:
   for line in name_sorted:
       fp.write(line+"\n")
