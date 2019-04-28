using KDTrees
#Usage julia Oversample.jl TA 20
Type = ARGS[1]
############################### Constantes
R_Smearing=int(ARGS[2])
D_Smearing=2*sin(deg2rad(R_Smearing/2.))
const granularity=10
#const Type="Sim_TA"
#const Type="Sim_Auger"

#const Type="TA"
#const Type="Auger"

############################### Constantes
if Type == "Sim_TA" || Type == "Sim_Auger"
    Datos=readdlm("100K_evts.txt")
    RA_Datos=Datos[:,1]
    DEC_Datos=90-Datos[:,2]
    #E_Datos=Datos[:,3]
    N_Datos=length(RA_Datos)
end
if Type == "TA" || Type == "Auger"
    RA_Datos=readdlm("RA_Datos.dat")
    DEC_Datos=90-readdlm("DEC_Datos.dat")
    #E_Datos=readdlm("E_Datos.dat")
    N_Datos=length(RA_Datos)
end
const N_lon = 360
const N_lat = 180
Overlap_Datos=zeros(N_lat*granularity,N_lon*granularity)
# data grid
x_d=sin(deg2rad(DEC_Datos)).*cos(deg2rad(RA_Datos))
y_d=sin(deg2rad(DEC_Datos)).*sin(deg2rad(RA_Datos))
z_d=cos(deg2rad(DEC_Datos))
data=[transpose(x_d),transpose(y_d),transpose(z_d)]
tree = KDTree(data,reorder = false)
if Type == "TA" || Type=="Sim_TA"
    theta_max_deg=55
    Lat_Obs=39.2969 # Latitud del TA
    lat_min=Lat_Obs-theta_max_deg
end
if  Type == "Auger" || Type=="Sim_Auger"
    theta_max_deg=80
    Lat_Obs=-35.2068 # Latitud del Auger
    lat_max=Lat_Obs+theta_max_deg
end
# loop over space grid
for i = 1:N_lat*granularity  
    if (i%granularity == 0 )
        Declination=90.-i/granularity
        print(" Declination=$Declination")
    end
    Declination=90.-i/granularity
    if ((Type=="TA" || Type=="Sim_TA")   && Declination >= lat_min) 
        
        for j = 1:N_lon*granularity
            x=sin(deg2rad((i-.5)/granularity))*cos(deg2rad((j-.5)/granularity))
            y=sin(deg2rad((i-.5)/granularity))*sin(deg2rad((j-.5)/granularity))
            z=cos(deg2rad((i-.5)/granularity))
            a=inball(tree, [x,y,z], D_Smearing, true)
            Overlap_Datos[i,j]=length(a)
        end
    end
 #   if ((Type=="Auger" || Type=="Sim_Auger")  && Declination <= lat_max) 
    if ((Type=="Auger" || Type=="Sim_Auger")) 

        for j = 1:N_lon*granularity
            x=sin(deg2rad((i-.5)/granularity))*cos(deg2rad((j-.5)/granularity))
            y=sin(deg2rad((i-.5)/granularity))*sin(deg2rad((j-.5)/granularity))
            z=cos(deg2rad((i-.5)/granularity))
            a=inball(tree, [x,y,z], D_Smearing, true)
            Overlap_Datos[i,j]=length(a)
        end
    end
end

print(" Granularity= $granularity N_Datos=$N_Datos  N_lon=$N_lon N_lat=$N_lat \n")
if Type == "Sim_TA" || Type == "Sim_Auger"
    writedlm("N_off.dat", Overlap_Datos)       
end
if (Type == "TA" || Type == "Auger") 
    writedlm("N_on.dat", Overlap_Datos)
    Val_Max = maximum(Overlap_Datos)
    print("Valor Maximo=$Val_Max \n") 
    Mat_Max=findn(Overlap_Datos .== Val_Max )
end
