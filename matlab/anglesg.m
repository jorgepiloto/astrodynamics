% ------------------------------------------------------------------------------
%
%                           function anglesg
%
%  this function solves the problem of orbit determination using three
%    optical sightings.  the solution function uses the gaussian technique.
%    there are lots of debug statements in here to test various options.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  23 dec 2003
%   8 oct 2007
%
%  inputs          description                    range / units
%    re           - radius earth, sun, etc        km
%    mu           - grav param earth, sun etc     km3/s2
%    tu           - time unit                     sec
%    rtasc1       - right ascension #1            rad
%    rtasc2       - right ascension #2            rad
%    rtasc3       - right ascension #3            rad
%    decl1       - declination #1                rad
%    decl2       - declination #2                rad
%    decl3       - declination #3                rad
%    jd1          - julian date of 1st sighting   days from 4713 bc
%    jd2          - julian date of 2nd sighting   days from 4713 bc
%    jd3          - julian date of 3rd sighting   days from 4713 bc
%    rs           - ijk site position vector      km
%
%  outputs        :
%    r            - ijk position vector at t2     km
%    v            - ijk velocity vector at t2     km / s
%
%  locals         :
%    l1           - line of sight vector for 1st
%    l2           - line of sight vector for 2nd
%    l3           - line of sight vector for 3rd
%    tau          - taylor expansion series about
%                   tau ( t - to )
%    tausqr       - tau squared
%    t21t23       - (t2-t1) * (t2-t3)
%    t31t32       - (t3-t1) * (t3-t2)
%    i            - index
%    d            -
%    rho          - range from site to sat at t2  km
%    rhodot       -
%    dmat         -
%    rs1          - site vectors
%    rs2          -
%    rs3          -
%    earthrate    - velocity of earth rotation
%    p            -
%    q            -
%    oldr         -
%    oldv         -
%    f1           - f coefficient
%    g1           -
%    f3           -
%    g3           -
%    l2dotrs      -
%
%  coupling       :
%    mag          - magnitude of a vector
%    detrminant   - evaluate the determinant of a matrix
%    factor       - find roots of a polynomial
%    matmult      - multiply two matrices together
%    gibbs        - gibbs method of orbit determination
%    hgibbs       - herrick gibbs method of orbit determination
%    angl         - angl between two vectors
%
%  references     :
%    vallado       2007, 429-439
%
% [r2,v2] = anglesg ( decl1,decl2,decl3,rtasc1,rtasc2,rtasc3,jd1,jd2,jd3,rs1,rs2,rs3, re, mu );
% ------------------------------------------------------------------------------

function [r2, v2] = anglesg ( decl1,decl2,decl3,rtasc1,rtasc2, ...
                    rtasc3,jd1,jd2,jd3, rs1, rs2, rs3, re, mu, tu );

    % -------------------------  implementation   -------------------------
    ddpsi = 0.0;  % delta correctinos not needed for this level of precision
    ddeps = 0.0;
    magr1in = 2.0*re; % initial guesses
    magr2in = 2.01*re;
    direct = 'y';  % direction of motion short way

    % ---------- set middle to 0, find decls to others -----------
    tau1= (jd1-jd2)*tu;  % sec
    tau3= (jd3-jd2)*tu;
    fprintf(1,'jd123 %14.6f %14.6f %14.6f tau %11.7f  %11.7f  \n',jd1,jd2,jd3,tau1,tau3);

    % ----------------  find line of sight unit vectors  ---------------
    l1(1)= cos(decl1)*cos(rtasc1);
    l1(2)= cos(decl1)*sin(rtasc1);
    l1(3)= sin(decl1);

    l2(1)= cos(decl2)*cos(rtasc2);
    l2(2)= cos(decl2)*sin(rtasc2);
    l2(3)= sin(decl2);

    l3(1)= cos(decl3)*cos(rtasc3);
    l3(2)= cos(decl3)*sin(rtasc3);
    l3(3)= sin(decl3);

    % ------------- find l matrix and determinant -----------------
    l1
    vs = [0 0 0]';
    aecef = [0 0 0]';
    %[l1eci,vs3,aeci] = ecef2eci(l1',vs,aecef,(jd1-2451545.0)/36525.0,jd1,0.0,0.0,0.0,0,ddpsi,ddeps);
    %[l2eci,vs3,aeci] = ecef2eci(l2',vs,aecef,(jd2-2451545.0)/36525.0,jd2,0.0,0.0,0.0,0,ddpsi,ddeps);
    %[l3eci,vs3,aeci] = ecef2eci(l3',vs,aecef,(jd3-2451545.0)/36525.0,jd3,0.0,0.0,0.0,0,ddpsi,ddeps);

    l1eci = l1;
    l2eci = l2;
    l3eci = l3;
    % leave these as they come since the topoc radec are alrady eci
    l1eci
    % --------- called lmati since it is only used for determ -----
    for i= 1 : 3
        lmatii(i,1) = l1eci(i);
        lmatii(i,2) = l2eci(i);
        lmatii(i,3) = l3eci(i);
        rsmat(i,1)  = rs1(i);
        rsmat(i,2)  = rs2(i);
        rsmat(i,3)  = rs3(i);
    end;
    lmatii
    
    fprintf(1,'rsmat eci %11.7f  %11.7f  %11.7f km \n',rsmat');

    % the order is right, but to print out, need '
    fprintf(1,'rsmat eci %11.7f  %11.7f  %11.7f \n',rsmat'/re);

    lmatii
    fprintf(1,'this should be the inverse of what the code finds later\n');
    li = inv(lmatii)
    d= det(lmatii);
    d
    % ------------------ now assign the inverse -------------------
    lmati(1,1) = ( l2eci(2)*l3eci(3)-l2eci(3)*l3eci(2)) / d;
    lmati(2,1) = (-l1eci(2)*l3eci(3)+l1eci(3)*l3eci(2)) / d;
    lmati(3,1) = ( l1eci(2)*l2eci(3)-l1eci(3)*l2eci(2)) / d;
    lmati(1,2) = (-l2eci(1)*l3eci(3)+l2eci(3)*l3eci(1)) / d;
    lmati(2,2) = ( l1eci(1)*l3eci(3)-l1eci(3)*l3eci(1)) / d;
    lmati(3,2) = (-l1eci(1)*l2eci(3)+l1eci(3)*l2eci(1)) / d;
    lmati(1,3) = ( l2eci(1)*l3eci(2)-l2eci(2)*l3eci(1)) / d;
    lmati(2,3) = (-l1eci(1)*l3eci(2)+l1eci(2)*l3eci(1)) / d;
    lmati(3,3) = ( l1eci(1)*l2eci(2)-l1eci(2)*l2eci(1)) / d;
    lmati
    lir = lmati*rsmat;

    % ------------ find f and g series at 1st and 3rd obs ---------
    %   speed by assuming circ sat vel for udot here ??
    %   some similartities in 1/6t3t1 ...
    % --- keep separated this time ----
    a1 =  tau3 / (tau3 - tau1);
    a1u=  (tau3*((tau3-tau1)^2 - tau3*tau3 )) / (6.0*(tau3 - tau1));
    a3 = -tau1 / (tau3 - tau1);
    a3u= -(tau1*((tau3-tau1)^2 - tau1*tau1 )) / (6.0*(tau3 - tau1));

    fprintf(1,'a1/a3 %11.7f  %11.7f  %11.7f  %11.7f \n',a1,a1u,a3,a3u );

    % --- form initial guess of r2 ----
    dl1=  lir(2,1)*a1 - lir(2,2) + lir(2,3)*a3;
    dl2=  lir(2,1)*a1u + lir(2,3)*a3u;
    dl1
    dl2

    % ------- solve eighth order poly not same as laplace ---------
    magrs2 = mag(rs2);
    l2dotrs= dot( l2,rs2 );
    fprintf(1,'magrs2 %11.7f  %11.7f  \n',magrs2,l2dotrs );

    poly( 1)=  1.0;  % r2^8th variable%%%%%%%%%%%%%%
    poly( 2)=  0.0;
    poly( 3)=  -(dl1*dl1 + 2.0*dl1*l2dotrs + magrs2^2);
    poly( 4)=  0.0;
    poly( 5)=  0.0;
    poly( 6)=  -2.0*mu*(l2dotrs*dl2 + dl1*dl2);
    poly( 7)=  0.0;
    poly( 8)=  0.0;
    poly( 9)=  -mu*mu*dl2*dl2;
    fprintf(1,'%11.7f \n',poly);
    
    rootarr = roots( poly );
    rootarr
    %fprintf(1,'rootarr %11.7f \n',rootarr);

    % ------------------ select the correct root ------------------
    bigr2= -99999990.0;
    % change from 1
    for j= 1 : 8
        if ( rootarr(j) > bigr2 ) & ( isreal(rootarr(j)) )
            bigr2= rootarr(j);
        end  % if (
    end
    bigr2

    % ------------ solve matrix with u2 better known --------------
    u= mu / ( bigr2*bigr2*bigr2 );

    c1= a1 + a1u*u;
    c2 = -1.0;
    c3= a3 + a3u*u;

    fprintf(1,'u %17.14f c1 %11.7f c3 %11.7f %11.7f \n',u,c1,c2,c3);

    cmat(1,1)= -c1;
    cmat(2,1)= -c2;
    cmat(3,1)= -c3;
    rhomat = lir*cmat;

    rhoold1=  rhomat(1,1)/c1;
    rhoold2=  rhomat(2,1)/c2;
    rhoold3=  rhomat(3,1)/c3;
    fprintf(1,'rhoold %11.7f %11.7f %11.7f \n',rhoold1,rhoold2,rhoold3);
    %   fprintf(1,'rhoold %11.7f %11.7f %11.7f \n',rhoold1/re,rhoold2/re,rhoold3/re);

    for i= 1 : 3
        r1(i)=  rhomat(1,1)*l1eci(i)/c1 + rs1(i);
        r2(i)=  rhomat(2,1)*l2eci(i)/c2 + rs2(i);
        r3(i)=  rhomat(3,1)*l3eci(i)/c3 + rs3(i);
    end
    fprintf(1,'r1 %11.7f %11.7f %11.7f \n',r1);
    fprintf(1,'r2 %11.7f %11.7f %11.7f \n',r2);
    fprintf(1,'r3 %11.7f %11.7f %11.7f \n',r3);

    pause;
    % -------- loop through the refining process ------------  while () for
    fprintf(1,'now refine the answer \n');
    rho2 = 999999.9;
    ll   = 0;
    while ((abs(rhoold2-rho2) > 1.0e-12) && (ll <= 0 ))  % ll <= 15
        ll = ll + 1;
        fprintf(1, ' iteration #%3i \n',ll );
        rho2 = rhoold2;  % reset now that inside while loop

        % ---------- now form the three position vectors ----------
        for i= 1 : 3
            r1(i)=  rhomat(1,1)*l1eci(i)/c1 + rs1(i);
            r2(i)= -rhomat(2,1)*l2eci(i)    + rs2(i);
            r3(i)=  rhomat(3,1)*l3eci(i)/c3 + rs3(i);
        end
        magr1 = mag( r1 );
        magr2 = mag( r2 );
        magr3 = mag( r3 );

        [v2,theta,theta1,copa,error] = gibbsh(r1,r2,r3, re, mu);

        rad = 180.0/pi;
        fprintf(1,'r1 %11.7f %11.7f %11.7f %11.7f %11.7f \n',r1,theta*rad,theta1*rad);
        fprintf(1,'r2 %11.7f %11.7f %11.7f \n',r2);
        fprintf(1,'r3 %11.7f %11.7f %11.7f \n',r3);
        fprintf(1,'w gibbs km/s       v2 %11.7f %11.7f %11.7f \n',v2);

        if ( (strcmp(error, '          ok') == 0) && (copa < 1.0/rad) ) % 0 is false
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2,v2, re, mu);
            fprintf(1,'coes init ans %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n',...
                p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
            % --- hgibbs to get middle vector ----
            [v2,theta,theta1,copa,error] = hgibbs(r1,r2,r3,jd1,jd2,jd3);
            fprintf(1,'using hgibbs: ' );
        end

        [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2,v2, re, mu);
        fprintf(1,'coes init ans %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n',...
            p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
        %fprintf(1,'dr %11.7f m %11.7f m/s \n',1000*mag(r2-r2ans),1000*mag(v2-v2ans) );

        if ( ll <= 8 )  % 4
            % --- now get an improved estimate of the f and g series --
            u= mu / ( magr2*magr2*magr2 );
            rdot= dot(r2,v2)/magr2;
            udot= (-3.0*mu*rdot) / (magr2^4);

            fprintf(1,'u %17.15f rdot  %11.7f udot %11.7f \n',u,rdot,udot );
            tausqr= tau1*tau1;
            f1=  1.0 - 0.5*u*tausqr -(1.0/6.0)*udot*tausqr*tau1;
            %                       - (1.0/24.0) * u*u*tausqr*tausqr
            %                       - (1.0/30.0)*u*udot*tausqr*tausqr*tau1;
            g1= tau1 - (1.0/6.0)*u*tau1*tausqr - (1.0/12.0) * udot*tausqr*tausqr;
            %                       - (1.0/120.0)*u*u*tausqr*tausqr*tau1
            %                       - (1.0/120.0)*u*udot*tausqr*tausqr*tausqr;
            tausqr= tau3*tau3;
            f3=  1.0 - 0.5*u*tausqr -(1.0/6.0)*udot*tausqr*tau3;
            %                       - (1.0/24.0) * u*u*tausqr*tausqr
            %                       - (1.0/30.0)*u*udot*tausqr*tausqr*tau3;
            g3= tau3 - (1.0/6.0)*u*tau3*tausqr - (1.0/12.0) * udot*tausqr*tausqr;
            %                       - (1.0/120.0)*u*u*tausqr*tausqr*tau3
            %                       - (1.0/120.0)*u*udot*tausqr*tausqr*tausqr;
            fprintf(1,'f1 %11.7f g1 %11.7f f3 %11.7f g3 %11.7f \n',f1,g1,f3,g3 );
        else
            % -------- use exact method to find f and g -----------
            theta  = angl( r1,r2 );
            theta1 = angl( r2,r3 );

            f1= 1.0 - ( (magr1*(1.0 - cos(theta)) / p ) );
            g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );  % - angl because backwards
            f3= 1.0 - ( (magr3*(1.0 - cos(theta1)) / p ) );
            g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );

        end

        c1=  g3 / (f1*g3 - f3*g1);
        c3= -g1 / (f1*g3 - f3*g1);

        fprintf(1,' c1 %11.7f c3 %11.7f %11.7f \n',c1,c2,c3);

        % ----- solve for all three ranges via matrix equation ----
        cmat(1,1)= -c1;
        cmat(2,1)= -c2;
        cmat(3,1)= -c3;
        rhomat = lir*cmat;

        fprintf(1,'rhomat %11.7f %11.7f %11.7f \n',rhomat);
        %        fprintf(1,'rhomat %11.7f %11.7f %11.7f \n',rhomat/re);

        rhoold1=  rhomat(1,1)/c1;
        rhoold2=  rhomat(2,1)/c2;
        rhoold3=  rhomat(3,1)/c3;
        fprintf(1,'rhoold %11.7f %11.7f %11.7f \n',rhoold1,rhoold2,rhoold3);
        %   fprintf(1,'rhoold %11.7f %11.7f %11.7f \n',rhoold1/re,rhoold2/re,rhoold3/re);

        for i= 1 : 3
            r1(i)=  rhomat(1,1)*l1eci(i)/c1 + rs1(i);
            r2(i)=  rhomat(2,1)*l2eci(i)/c2 + rs2(i);
            r3(i)=  rhomat(3,1)*l3eci(i)/c3 + rs3(i);
        end
        fprintf(1,'r1 %11.7f %11.7f %11.7f \n',r1);
        fprintf(1,'r2 %11.7f %11.7f %11.7f \n',r2);
        fprintf(1,'r3 %11.7f %11.7f %11.7f \n',r3);

        fprintf(1,'====================next loop \n');
        % ----------------- check for convergence -----------------
        pause
        fprintf(1,'rhoold while  %16.14f %16.14f \n',rhoold2,rho2);
    end   % while the ranges are still changing

    % ---------------- find all three vectors ri ------------------
    for i= 1 : 3
        r1(i)=  rhomat(1,1)*l1eci(i)/c1 + rs1(i);
        r2(i)= -rhomat(2,1)*l2eci(i)    + rs2(i);
        r3(i)=  rhomat(3,1)*l3eci(i)/c3 + rs3(i);
    end

