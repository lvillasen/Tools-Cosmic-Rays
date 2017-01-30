import pandas as pd
import numpy as np
import cosmolopy.distance as cd
import astroconvert
cosmo = {'omega_M_0':0.286, 'omega_lambda_0':0.714, 'omega_k_0':0.0, 'h':0.696}
VCV=pd.read_csv('BAT_70m_catalog_20nov2012.txt', sep='|', skiprows=2)
#distance=np.zeros(len(VCV));RA_AGN=np.zeros(len(VCV));DEC_AGN=np.zeros(len(VCV))
distance=[];RA_AGN=[];DEC_AGN=[];name=[];flux=[];redshift=[]
j=0
for i in range(len(VCV)):
    if VCV.values[i][17] >0 and VCV.values[i][17] <.0235:
        RA2=float("{0:.2f}".format(float(VCV.values[i][2])))
        DEC2=float("{0:.2f}".format(float(VCV.values[i][3])))
        RA_AGN.append(RA2)
        DEC_AGN.append(DEC2)
        redshift.append(VCV.values[i][17])
        D=cd.comoving_distance(VCV.values[i][17], **cosmo)
        D2=float("{0:.2f}".format(float(D)))
        distance.append(D2)
        name.append(VCV.values[i][5])
        flux.append(VCV.values[i][9])
        j+=1
sorted_index, sorted_d = zip(*sorted([(i,e) for i,e in enumerate(flux)], key=lambda x:x[1]))
sorted_index=sorted_index[::-1]
RA_AGN_sorted=[RA_AGN[i] for i in sorted_index]
DEC_AGN_sorted=[DEC_AGN[i] for i in sorted_index]
distance_sorted=[distance[i] for i in sorted_index]   
name_sorted=[name[i] for i in sorted_index] 
flux_sorted=[flux[i] for i in sorted_index] 
redshift_sorted=[redshift[i] for i in sorted_index]    


print  "Total number of objects=",len(VCV)
print  "Total number of objects with measured redshift < 0.0235=",j

print "Name","\t","RA","\t","DEC","\t","\t","Flux","\t","Redshift"
for i in range(20):
    print name_sorted[i],"\t", RA_AGN_sorted[i],"\t", DEC_AGN_sorted[i],"\t",flux_sorted[i],"\t",redshift_sorted[i]
print " ......."
print name_sorted[-1],"\t", RA_AGN_sorted[-1],"\t", DEC_AGN_sorted[-1],"\t",flux_sorted[-1],"\t",redshift_sorted[-1]
np.savetxt("RA_AGN.dat", RA_AGN_sorted)
np.savetxt("DEC_AGN.dat", DEC_AGN_sorted)
np.savetxt("distance.dat", distance_sorted)
np.savetxt("flux.dat", flux_sorted)
with open ("name.dat","w") as f:
   for line in name_sorted:
       f.write(line+"\n")
