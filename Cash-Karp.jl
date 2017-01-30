
# Cash-Karp parameters
A = [ 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 ]
B=zeros(6,6)
B[2,:]=[0.2,0,0,0,0,0]
B[3,:]=[3.0/40.0, 9.0/40.0,0,0,0,0]
B[4,:]=[0.3, -0.9, 1.2,0,0,0]
B[5,:]=[-11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0,0,0]
B[6,:]=[1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0,0]
C  = [37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0]
DC = [C[1]-2825.0/27648.0, C[2]-0.0, C[3]-18575.0/48384.0,C[4]-13525.0/55296.0, C[5]-277.00/14336.0, C[6]-0.25]
function Cash_Karp(f,XX0,d_max)
    c=2.998*10^(8.); # Usamos MKS
    kpc=3.0857*10^19.;
    kyear=3.154*10^10.;
    t=0;;niter=0;niter_si_avanza=0;niter_no_avanza=0
    
    hmax=0.50*kpc/c # 400 pc/c
    delta_error_max=.0001 
    D=0;
    K=zeros(6,6)
    h=hmax
    n_trial=0
    while (D<d_max*kpc )
        if D<=20*kpc && n_trial==0;
            h=hmax
        elseif D>20*kpc && n_trial==0;
            h=hmax*1000
        end
        D=sqrt(XX0[1]*XX0[1]+XX0[2]*XX0[2]+XX0[3]*XX0[3]);
        K[:,1]=h*f(XX0,t)
        K[:,2]=h*f(XX0+B[2,1]*K[:,1],t+h*A[2])
        K[:,3]=h*f(XX0+B[3,2]*K[:,2]+B[3,1]*K[:,1],t+h*A[3])
        K[:,4]=h*f(XX0+B[4,3]*K[:,3]+B[4,2]*K[:,2]+B[4,1]*K[:,1],t+h*A[4])
        K[:,5]=h*f(XX0+B[5,4]*K[:,4]+B[5,3]*K[:,3]+B[5,2]*K[:,2]+B[5,1]*K[:,1],t+h*A[5])
        K[:,6]=h*f(XX0+B[6,5]*K[:,5]+B[6,4]*K[:,4]+B[6,3]*K[:,3]+B[6,2]*K[:,2]+B[6,1]*K[:,1],t+h*A[6])
        
        desp=maximum(C[1]*K[1:3,1]+ C[2]*K[1:3,2]+ C[3]*K[1:3,3]+ C[4]*K[1:3,4]+ C[5]*K[1:3,5]+ C[6]*K[1:3,6])
        error=maximum(DC[1]*K[1:3,1]+ DC[2]*K[1:3,2]+ DC[3]*K[1:3,3]+ DC[4]*K[1:3,4]+ DC[5]*K[1:3,5]+ DC[6]*K[1:3,6])
        delta_error=abs(error/desp)
        niter += 1;
        if delta_error < delta_error_max
            XX0=XX0 + C[1]*K[:,1]+ C[2]*K[:,2]+ C[3]*K[:,3]+ C[4]*K[:,4]+ C[5]*K[:,5]+ C[6]*K[:,6]
            t += h;niter_si_avanza+=1
            h=hmax;n_trial=0
        else
            h=h/2.; niter_no_avanza+=1;n_trial=1
        end
    end 
    return XX0,t,niter,niter_si_avanza 
end