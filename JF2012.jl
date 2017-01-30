#spiral arm parameters
function logisticFunction(x,x0, w)
    return 1. / (1. + exp(-2. * (abs(x) - x0) / w));
end



function getRegularField(x,y,z)
    kpc=1.0;muG=1.0; 
    # spiral arm parameters
    pitch = 11.5 * pi*1.0 / 180;
    sinPitch = sin(pitch);
    cosPitch = cos(pitch);
    tan90MinusPitch = tan(pi / 2 - pitch);
    rArms=zeros(8)
    rArms[1] = 5.1 * kpc;
    rArms[2] = 6.3 * kpc;
    rArms[3] = 7.1 * kpc;
    rArms[4] = 8.3 * kpc;
    rArms[5] = 9.8 * kpc;
    rArms[6] = 11.4 * kpc;
    rArms[7] = 12.7 * kpc;
    rArms[8] = 15.5 * kpc;

    # regular field parameters
    bRing = 0.1 * muG;
    hDisk = 0.40 * kpc;
    wDisk = 0.27 * kpc;
    bDisk=zeros(8)
    bDisk[1] = 0.1 * muG;
    bDisk[2] = 3.0 * muG;
    bDisk[3] = -0.9 * muG;
    bDisk[4] = -0.8 * muG;
    bDisk[5] = -2.0 * muG;
    bDisk[6] = -4.2 * muG;
    bDisk[7] = 0.0 * muG;
    bDisk[8] = 2.7 * muG;


    bNorth = 1.4 * muG;
    bSouth = -1.1 * muG;
    rNorth = 9.22 * kpc;
    rSouth = 17 * kpc;
    wHalo = 0.20 * kpc;
    z0 = 5.3 * kpc;

    bX = 4.6 * muG;
    thetaX0 = 49.0 * pi*1.0 / 180;
    sinThetaX0 = sin(thetaX0);
    cosThetaX0 = cos(thetaX0);
    tanThetaX0 = tan(thetaX0);
    rXc = 4.8 * kpc;
    rX = 2.9 * kpc;

    

    # turbulent field parameters
    bDiskTurb=zeros(8)
    bDiskTurb[1] = 10.81 * muG;
    bDiskTurb[2] = 6.96 * muG;
    bDiskTurb[3] = 9.59 * muG;
    bDiskTurb[4] = 6.96 * muG;
    bDiskTurb[5] = 1.96 * muG;
    bDiskTurb[6] = 16.34 * muG;
    bDiskTurb[7] = 37.29 * muG;
    bDiskTurb[8] = 10.35 * muG;

    bDiskTurb5 = 7.63 * muG;
    zDiskTurb = 0.61 * kpc;

    bHaloTurb = 4.68 * muG;
    rHaloTurb = 10.97 * kpc;
    zHaloTurb = 2.84 * kpc;
    bx=0;by=0;bz=0;
    r = sqrt(x * x + y * y); # in-plane radius
    d = sqrt(x * x + y * y + z*z); # distance to galactic center
    if (d < 1 ) || (d > 20 )
        return [bx,by,bz]; # 0 field for d < 1 kpc or d > 20 kpc
    end

    phi=atan2(y,x); # phi de -pi a pi
    sinPhi = sin(phi);
    cosPhi = cos(phi);

    lfDisk = logisticFunction(z, hDisk, wDisk);

    # disk field
    if r > 3  
        if r < 5  
            # molecular ring
            bMag = bRing * (5 * kpc / r) * (1 - lfDisk);
            bx += -bMag * sinPhi;
            by += bMag * cosPhi;
        else 
            # spiral region
            i_0=7            
            
            for i=1:7
                r11=rArms[i] *exp((phi*pitch-pi*pitch)) ;
                r12=rArms[i+1] *exp((phi*pitch-pi*pitch)) ;
                r21=rArms[i] *exp(((phi+2*pi)*pitch-pi*pitch));
                r22=rArms[i+1] *exp(((phi+2*pi)*pitch-pi*pitch));
                r31=rArms[i] *exp(((phi-2*pi)*pitch-pi*pitch));
                r32=rArms[i+1] *exp(((phi-2*pi)*pitch-pi*pitch));
                if (r>=r11 && r<r12) ||  (r>=r21 && r<r22) ||  (r>=r31 && r<r32) 
                    i_0=i
                end
            end
            bMag = bDisk[i_0];

            bMag *= (5 * kpc / r) * (1 - lfDisk);
            bx += bMag * (sinPitch * cos(phi) - cosPitch * sin(phi));
            by += bMag * (sinPitch * sin(phi) + cosPitch * cos(phi));  
        end
    end
            
        
# toroidal halo field
    bMagH = exp(-abs(z) / z0) * lfDisk;
    if (z >= 0)
        bMagH *= bNorth * (1 - logisticFunction(r, rNorth, wHalo));
    else
        bMagH *= bSouth * (1 - logisticFunction(r, rSouth, wHalo));
    end
    bx += -bMagH * sinPhi;
    by += bMagH * cosPhi;
    
    
# poloidal halo field

    rc = rXc + abs(z) / tanThetaX0;
    if (r < rc) 
        # varying elevation region
        rp = r * rXc / rc;
        bMagX = bX * exp(-1 * rp / rX) * (rp / r)^ 2.;
        thetaX = atan2(abs(z), (r - rp));
        if (z == 0) 
            thetaX = pi / 2.;
        end
        sinThetaX = sin(thetaX);
        cosThetaX = cos(thetaX);
    else 
        # constant elevation region
        rp = r - abs(z) / tanThetaX0;
        bMagX = bX * exp(-rp / rX) * (rp / r);
        sinThetaX = sinThetaX0;
        cosThetaX = cosThetaX0;
    end
    
    bx += sign(z) * bMagX * cosThetaX * cosPhi;
    by += sign(z) * bMagX * cosThetaX * sinPhi;
    bz += bMagX * sinThetaX;
    #print ("JF2012 x=$x y=$y z=$z \n");
    #print ("JF2012 Bx=$bx By=$by Bz=$bz \n");
    return [bx,by,bz];  
end
function setStriatedGrid()
    N = 400;
    StriatedGrid=zeros(N,N,N)
    for ix = 1:N
        for iy = 1:N
            for iz = 1:N
                StriatedGrid[ix,iy,iz]= 2*round(Int,rand())-1
            end
        end
    end
    return StriatedGrid
end
function getStriatedField(x,y,z,StriatedGrid)
    # striated field parameter
    sqrtbeta = sqrt(1.36);
    N=400
    inside=1.0
    D=sqrt(x^2+y^2+z^2)
    if (D>10.0); inside=0.;end
    ix_grid=round(Int,10*(x+20.0)+.5)
    #if (ix_grid < 1 || ix_grid > N) ;print("ix_grid=",ix_grid) ;end
    if ix_grid < 1;ix_grid=1 ;end
    if ix_grid > N;ix_grid=N ;end
    iy_grid=round(Int,10*(y+20.0)+.5)
    #if (iy_grid < 1 || iy_grid > N) ;print("iy_grid=",iy_grid) ;end
    if iy_grid < 1 ; ;iy_grid=1 ;end
    if iy_grid > N;iy_grid=N ;end
    iz_grid=round(Int,10*(z+20.0)+.5)
    #if (iz_grid < 1 || iz_grid > N) ;print("iz_grid=",iz_grid) ;end
    if iz_grid < 1; ;iz_grid=1 ;end
    if iz_grid > N;iz_grid=N ;end
    (bx_S,by_S,bz_S)=inside*sqrtbeta*StriatedGrid[ix_grid,iy_grid,iz_grid]*getRegularField(x,y,z);
    return [bx_S,by_S,bz_S];  
end
