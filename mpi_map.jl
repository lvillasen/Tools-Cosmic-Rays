include("JF2012.jl")
include("Cash-Karp.jl")
include("Propa_JF2012.jl")

using MPI
using SkyCoords
#Pkg.add("SkyCoords")


MPI.Init()
comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const size = MPI.Comm_size(comm)
for k = 1:size
    if  k-1==rank
        print("Hello World, I am process $rank of $size running on name. \n" )
    end
end
#const name = MPI.Comm_name(comm)
############################### Constantes
const Coord_EQ = 1 # 1 for map in EQUATORIAL coord.; 0 for GALACTIC coord.
const granularity=5
const Charge=26 # 1 for p, 8 for O, 26 for Fe
const Mass =56 # 1 for p, 16 for O, 56 for Fe
const Energy =60
############################### Constants
const N_lon = 360
const N_lat = 180
Desviacion=zeros(N_lat*granularity,N_lon*granularity)
 
   
for i = 1:N_lat*granularity
    if  (i-1)%size == rank
        print (" i=$i ")
        for j = 1:N_lon*granularity
            if Coord_EQ == 1
                eq_dec=90.-(i-.5)/granularity
                eq_ra= (j-.5)/granularity
                c1 = FK5Coords{2000}(deg2rad(eq_ra), deg2rad(eq_dec))  # inputs are ra, dec in radians
                c2 = convert(GalCoords, c1)
                galactic_l=rad2deg(c2.l)
                galactic_b=rad2deg(c2.b)
                alpha=Propa_JF2012_simple(Charge,Mass,Energy, galactic_l, galactic_b,0,true,false,false,20000.);
                Desviacion[i,j]=alpha
            else
                galactic_b=90.-(i-.5)/granularity
                galactic_l= (j-.5)/granularity-180.
                alpha=Propa_JF2012_simple(Charge,Mass,Energy, galactic_l, galactic_b,0,true,false,false,20000.);
                Desviacion[i,j]=alpha
            end
        end
    
    end    
end
data2=MPI.Reduce(Desviacion, MPI.SUM, 0, comm)

if rank==0
    print ("Granularity= $granularity  N_lon=$N_lon N_lat=$N_lat \n")
    data21 = reshape (data2, N_lat*granularity,N_lon*granularity)
    if Coord_EQ == 1
        writedlm("Desviacion_EQ_Fe60g5.dat", data21) 
    else
        writedlm("Desviacion_GAL_p60g5.dat", data21) 

    end
end

MPI.Finalize()  
