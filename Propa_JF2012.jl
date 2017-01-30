
function Propa_JF2012(Charge,Mass,Energy, galactic_l, galactic_b,verbose,Regular,Striated,EG,d_max)
    # Entradas:  Charge,Mass,Energy, galactic_l, galactic_b, verbose
    

    ######## Constantes
    EeV=1.602*10^(-1.);
    Mp=1.67*10^(-27.);
    Qp=1.602*10^(-19.) 
    c=2.998*10^(8.);
    kpc=3.0857*10^19.;
    muGauss=10.^(-10.); # Tesla
    kyear=3.154*10^10.;

    ############## Cond. iniciales
    E0=Energy*EeV;# p=Ev/c^2
    q=-Charge*Qp; # units of p charge, negative foe backpropagation
    m0=Mass*Mp # units of p mass
    
    ############## Cond. iniciales
    galactic_l_rad=deg2rad(galactic_l);
    galactic_b_rad=deg2rad(galactic_b);
    nx=cos(galactic_b_rad)*cos(galactic_l_rad);
    ny=cos(galactic_b_rad)*sin(galactic_l_rad);
    nz=sin(galactic_b_rad);


    n_inicial=[nx,ny,nz]
    n_inicial=n_inicial/norm(n_inicial);
    v0=c*sqrt(1-m0^2c^4/E0^2)*n_inicial;
    x0=-8.5 * kpc; # Sun position
    y0=0;
    z0=0;
    #X0_noRel=[x0,y0,z0,v0[1],v0[2],v0[3]];
    X0_Rel=[x0,y0,z0,E0*v0[1]/c^2,E0*v0[2]/c^2,E0*v0[3]/c^2];
    if Striated;
        StriatedGrid=setStriatedGrid()
    end  
    

    function f_Rel(X,t)
        # X es [x,y,z,px,py,pz]
        # Distancia en m y tiempo en s
        c=2.998*10^(8.);
        M=m0*sqrt(1+(X[4]^2+X[5]^2+X[6]^2)/(m0*c)^2); #  masa relativista 
        #M=E0/c^2 # masa constante
        (Bx,By,Bz)=(0.,0.,0.)
        if Regular
            (Bx,By,Bz)=getRegularField(X[1]/kpc,X[2]/kpc,X[3]/kpc)*muGauss; # B en Teslas
        end
        if Striated; 
            B_S=getStriatedField(X[1]/kpc,X[2]/kpc,X[3]/kpc,StriatedGrid)*muGauss; # B en Teslas
            Bx+=B_S[1];By+=B_S[2];Bz+=B_S[3]
        end  
        if EG;
            B_sgal=EGMF(X[1]/kpc,X[2]/kpc,X[3]/kpc)*muGauss;# B en Teslas
            Bx+=B_sgal[1];Bx+=B_sgal[2];Bx+=B_sgal[3];
        end 
        
        return [X[4]/M,X[5]/M,X[6]/M,q/M*(X[5]*Bz-X[6]*By),q/M*(X[6]*Bx-X[4]*Bz),q/M*(X[4]*By-X[5]*Bx)]
    end


    #xt, t_f, niter=RG4_D(f_Rel,X0_Rel);
    xt, t_f, niter, niter_advance=Cash_Karp(f_Rel,X0_Rel,d_max);
    x_f=xt[1]/kpc;
    y_f=xt[2]/kpc;
    z_f=xt[3]/kpc;
    px=xt[4];
    py=xt[5];
    pz=xt[6];
    p_f=[px,py,pz];
    M=m0*sqrt(1+(px^2+py^2+pz^2)/(m0*c)^2);
    E_f=M*c^2/EeV;
    t_f=t_f/kyear;
    R_f=sqrt(x_f^2+y_f^2+z_f^2);
    alpha=rad2deg(acos(dot(n_inicial,p_f)/norm(p_f)));
    l_f=rad2deg(atan2(py,px));
    if l_f <0; l_f += 360.;end
    b_f=rad2deg(atan2(pz,sqrt(px*px+py*py)));
    if verbose ==1
        @printf("niter \t niter_advance \tE_f\tx_f\ty_f\tz_f\tl_f \tb_f \talpha \tt_f \n");
        @printf("%d \t %d \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",niter,niter_advance,E_f,x_f,y_f,z_f,l_f, b_f,alpha,t_f);  
    end
    alpha=@sprintf("%.2f",alpha)
    alpha=float(alpha)
    return alpha
end