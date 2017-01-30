const theta_min_deg=0.0
const theta_min=deg2rad(theta_min_deg)

#const Type="Auger"
const Type="TA"
const N_Datos=100000
const N_lat=180
const N_lon=360

CienMil_evts=zeros(N_Datos,3)
#####################
if Type=="TA"
    const Lat_Obs=39.2969
    const theta_max_deg=55.0
elseif Type=="Auger"
    const Lat_Obs=-35.2068
    const theta_max_deg=80.0
end
const theta_max=deg2rad(theta_max_deg)
const Theta_0=deg2rad(Lat_Obs-90.)
Ratio_vert_to_inclined=(sin(deg2rad(60))^2)/(sin(deg2rad(80))^2-sin(deg2rad(60))^2)
print("Ratio_vert_to_inclined=$Ratio_vert_to_inclined\n")
# PARA EL AUGER LA RAZON DE EXPOSURES DE VERT/INCLINADO = 3.52
Ratio_vert_to_inclined_data=3.52
# HAY QUE REDUCIR LOS EVENTOS INCLINADOS EN UN FACTOR ALFA TAL QUE R_solid_angle/ALFA=R_Exposure
# ES DECIR QUE ALFA=R_solid_angle/R_Exposure
alpha=Ratio_vert_to_inclined/Ratio_vert_to_inclined_data
print("alpha=$alpha")
N_Observed=1
while N_Observed<=N_Datos
    theta=rand(1)[1]*(theta_max-theta_min)+theta_min
        X_rand=rand(1)[1]
        theta_deg=rad2deg(theta) 
        Bool_1=((theta_deg>=60 && X_rand < alpha) || theta_deg<60)
        Bool_2=Bool_1 || (Tipo_de_Datos=="TA")
    if (rand(1)[1]<=sin(theta)*cos(theta)) && Bool_2
        phi=rand(1)[1]*2*pi-pi*1
        x1=sin(theta)*cos(phi)
        y1=sin(theta)*sin(phi)
        z1=cos(theta)
        #x1_p=x1*np.cos(Theta_0)-z1*np.sin(Theta_0)
        theta_p_calc=rad2deg(acos(sin(theta)*cos(phi)*sin(Theta_0) + cos(theta)*cos(Theta_0)))
        DEC=90-theta_p_calc
        RA=rand(1)[1]*N_lon
        CienMil_evts[N_Observed,1]= RA
        CienMil_evts[N_Observed,2]=DEC
        CienMil_evts[N_Observed,3]=theta_deg
        N_Observed+=1
    end
end
writedlm("100K_evts.txt",CienMil_evts)
