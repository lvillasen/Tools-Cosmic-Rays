# CODIGO PARA Calcular la significancia a la Li-Ma y P-value DE POISSON


using StatsFuns

############################### Constantes
const Type="TA"
#const Type="Auger"

const N_lon = 360
const N_lat = 180
const granularity=10
if Type=="TA"
    const eta=72./100000.
elseif Type=="Auger"
    const eta=201./100000.
end

N_on=readdlm("N_on.dat")
N_off=readdlm("N_off.dat")

S_LM=zeros(N_lat*granularity,N_lon*granularity)
P_VALUE=zeros(N_lat*granularity,N_lon*granularity)
n=0
for i = 1:N_lat*granularity   
    
    for j = 1:N_lon*granularity 
 
        if N_on[i,j] >0
            a=((1+eta)*N_on[i,j])/(eta*(N_on[i,j]+N_off[i,j]))
            val_1=N_on[i,j]*log(a)
        else
            val_1=0
        end
        if N_off[i,j] >0
            b=((1+eta)*N_off[i,j])/(N_on[i,j]+N_off[i,j])
            val_2=N_off[i,j]*log(b)
        else
            
            val_2=0
        end
        val=2*(val_1+val_2)
        if val<0
            n=n+1

            print("i=$i,j=$j,val=$val \n")

            
            val=0
        end
        c=sqrt(val)

        S_LM[i,j]=sign(N_on[i,j]-eta*N_off[i,j])*c
        
        
        if N_off[i,j] == 0.0
            P_VALUE[i,j]= 1
        else
            P_VALUE[i,j]=1-poiscdf(eta*N_off[i,j],N_on[i,j]-1)
        end
        
   end
end
writedlm("S_LM.dat", S_LM)
writedlm("P_VALUE.dat", -log10(P_VALUE))
#writedlm("S_LM_GAL.dat", S_LM)
Max_S=maximum(S_LM)
Max_P=maximum(-log10(P_VALUE))
@printf("Max Val S_LM=%.4f \n",Max_S)
@printf("Max Val -log10(P_VALUE)=%.4f \n",Max_P)
#print ("Min Val P_VALUE=$Min_P \n")
