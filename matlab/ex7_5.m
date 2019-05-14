%     -----------------------------------------------------------------
%
%                              Ex7_5.m
%
%  this file demonstrates example 7-5.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            30 mar 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
% calls
%  run4lamb  - runs the 4 cases and outputs results to a file
%            - each case has 250 minutes of a loop
%
%     *****************************************************************

    fid = 1;    
    directory = 'd:\codes\library\matlab\';
    outfile = fopen(strcat(directory,'tlambfig.out'), 'wt');
    
    constastro;

    % ---------------------------- book tests -----------------------------
    fprintf(1,'\n-------- lambert test book pg 497 short way \n' );
    r1 = [ 2.500000,    0.000000 ,   0.000000]*re;
    r2 = [ 1.9151111,   1.6069690,   0.000000]*re;
    % original orbit, assume circular
    v1 = [0 0 0];
    v1(2) = sqrt(mu/r1(1));
    ang = atan(r2(2)/r2(1));
    v2 = [-sqrt(mu/r2(2))*cos(ang);sqrt(mu/r2(1))*sin(ang);0.0];
    fprintf(1,'\n v1 %16.8f%16.8f%16.8f\n',v1 );
    fprintf(1,'\n v2 %16.8f%16.8f%16.8f\n',v2 );
   
    magr1 = mag(r1);
    magr2 = mag(r2);

    % this value stays constant in all calcs, vara changes with df
    cosdeltanu = dot(r1,r2) / (magr1*magr2);  
    
    dtsec = 76.0*60.0;
    dm = 's';
    df = 'd';
    nrev = 0;
    [ tbi, tbil] = lambgettbiu(r1, r2, 5);
    fprintf(1,' r1 %16.8f%16.8f%16.8f\n',r1 );
    fprintf(1,' r2 %16.8f%16.8f%16.8f\n',r2 );
    fprintf(1,'From universal variables \n%11.7f %11.7f s \n',tbi(1,1),tbi(1,2));
    fprintf(1,'%11.7f %11.7f s \n',tbi(2,1),tbi(2,2));
    fprintf(1,'%11.7f %11.7f s \n',tbi(3,1),tbi(3,2));
    fprintf(1,'%11.7f %11.7f s \n',tbi(4,1),tbi(4,2));
    fprintf(1,'%11.7f %11.7f s \n\n',tbi(5,1),tbi(5,2));
    fprintf(1,'%11.7f %11.7f s \n',tbil(1,1),tbil(1,2));
    fprintf(1,'%11.7f %11.7f s \n',tbil(2,1),tbil(2,2));
    fprintf(1,'%11.7f %11.7f s \n',tbil(3,1),tbil(3,2));
    fprintf(1,'%11.7f %11.7f s \n',tbil(4,1),tbil(4,2));
    fprintf(1,'%11.7f %11.7f s \n',tbil(5,1),tbil(5,2));

    [tbik, tbilk] = lambgettbik(r1, r2, 5);
    fprintf(1,'From k universal variables \n%11.7f %11.7f s %11.7f tu\n',tbik(1,1),tbik(1,2)*tusec,tbik(1,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n',tbik(2,1),tbik(2,2)*tusec,tbik(2,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n',tbik(3,1),tbik(3,2)*tusec,tbik(3,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n',tbik(4,1),tbik(4,2)*tusec,tbik(4,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n\n',tbik(5,1),tbik(5,2)*tusec,tbik(5,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n',tbilk(1,1),tbilk(1,2)*tusec,tbilk(1,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n',tbilk(2,1),tbilk(2,2)*tusec,tbilk(2,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n',tbilk(3,1),tbilk(3,2)*tusec,tbilk(3,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n',tbilk(4,1),tbilk(4,2)*tusec,tbilk(4,2));
    fprintf(1,'%11.7f %11.7f s %11.7f tu\n',tbilk(5,1),tbilk(5,2)*tusec,tbilk(5,2));

    [minenergyv, aminenergy, tminenergy, tminabs] = lambertmin ( r1, r2, 'd', 0 );
    fprintf(1,' minenergyv %16.8f %16.8f %16.8f a %11.7f  dt %11.7f  %11.7f \n', minenergyv, aminenergy, tminenergy, tminabs );
     
    [minenergyv, aminenergy, tminenergy, tminabs] = lambertmin ( r1, r2, 'r', 0 );
    fprintf(1,' minenergyv %16.8f %16.8f %16.8f a %11.7f  dt %11.7f  %11.7f \n', minenergyv, aminenergy, tminenergy, tminabs );

    
    dtwait = 0.0;
    fprintf(1,'\n-------- lambertu test \n' );
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, dm, df,nrev, dtwait, dtsec, tbi, fid );
    fprintf(1,' v1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    
    % run the 6 cases
       fprintf(1,' TEST ------------------ s  d  0 rev \n');
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, 's', 'd', 0, dtwait, 21300, tbi, fid );
    fprintf(1,' v1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    [v1t, v2t, errorl] = lambertk ( r1, v1, r2, 's', 'd', 0, dtwait, 21300, tbik, fid );
    fprintf(1,' v1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1dv, re, mu);
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'ans coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.3f \n',...
%             p,a,ecc,incl*rad,omega*rad,argp*rad, period );

       fprintf(1,' TEST ------------------ l  r  0 rev \n');
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, 'l', 'r', 0, dtwait, 21300, tbil, fid );
    fprintf(1,' v1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    [v1t, v2t, errorl] = lambertk ( r1, v1, r2, 'l', 'r', 0, dtwait, 21300, tbilk, fid );
    fprintf(1,' v1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
pause;
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1dv, re, mu);
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'ans coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.3f \n',...
%             p,a,ecc,incl*rad,omega*rad,argp*rad, period );

       fprintf(1,' TEST ------------------ s  d  1 rev \n');
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, 's', 'd', 1, dtwait, 21300, tbi, fid );
    fprintf(1,'uv1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    [v1t, v2t, errorl] = lambertk ( r1, v1, r2, 's', 'd', 1, dtwait, 21300, tbik, fid );
    fprintf(1,'kv1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1dv, re, mu);
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'ans coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.3f \n',...
%             p,a,ecc,incl*rad,omega*rad,argp*rad, period );
    
       fprintf(1,' TEST ------------------ s  r  1 rev \n');
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, 's', 'r', 1, dtwait, 21300, tbil, fid );
    fprintf(1,'uv1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    [v1t, v2t, errorl] = lambertk ( r1, v1, r2, 's', 'r', 1, dtwait, 21300, tbilk, fid );
    fprintf(1,'kv1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1dv, re, mu);
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'ans coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.3f \n',...
%             p,a,ecc,incl*rad,omega*rad,argp*rad, period );
    
       fprintf(1,' TEST ------------------ l  d  1 rev \n');
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, 'l', 'd', 1, dtwait, 21300, tbi, fid );
    fprintf(1,'uv1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    [v1t, v2t, errorl] = lambertk ( r1, v1, r2, 'l', 'd', 1, dtwait, 21300, tbik, fid );
    fprintf(1,'kv1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1dv, re, mu);
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'ans coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.3f \n',...
%             p,a,ecc,incl*rad,omega*rad,argp*rad, period );
    
       fprintf(1,' TEST ------------------ l  r  1 rev \n');
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, 'l', 'r', 1, dtwait, 21300, tbil, fid );
    fprintf(1,'uv1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    [v1t, v2t, errorl] = lambertk ( r1, v1, r2, 'l', 'r', 1, dtwait, 21300, tbilk, fid );
    fprintf(1,'kv1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1dv, re, mu);
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'ans coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.3f \n',...
%             p,a,ecc,incl*rad,omega*rad,argp*rad, period );

    fprintf(1,'\n-------- lambertb test \n' );
    [v1t, v2t, errorl] = lambertb ( r1, v1, r2, dm, nrev, dtsec );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2dv %16.8f%16.8f%16.8f\n',v2t );

    fprintf(1,'\n-------- lambertk test \n' );



    %    % create plot
%    run4lamb;
    % test Blair's cases
    r1 = [5690.2192300, 3309.6237700, 1311.3050400];
    r2 = [-5691.9147514,-3310.6095425, -1311.6957711];
    dtsec = 2750.0;
    p = 6713.05795271024;
    ecc = 0.00374086606300561;
    dnu = 179.999997168353/rad;
    
    [v1th, v2th] = lambhodograph( r1, v1, r2, p, ecc, dnu, dtsec );
    fprintf(1,'h v1t %16.8f %16.8f %16.8f %16.8f %16.8f \n',v1th, mag(v1th), dnu*rad );
    fprintf(1,'h v2t %16.8f %16.8f %16.8f %16.8f\n',v2th, mag(v2th) );
	fprintf(1,'answer    -3.67875419       6.71895632     -0.84722693 \n');
    [r3h,v3h] =  kepler  ( r1,v1th, dtsec );
    fprintf(1,'ANS hr3 %16.8f %16.8f %16.8f v3 %16.8f %16.8f %16.8f\n',r3h, v3h );
    
    fprintf(1,'------------------done with 180 transfers---------------------\n');
    pause;

    
    
   

nrev = 1;
lower = 4.0*nrev^2*pi*pi;
upper = 4.0*(nrev + 1.0)^2*pi*pi;
psiU = (lower+upper)/2*0.8;
[c2,c3] = findc2c3( psiU );
c2dotU = 1.0/(2*psiU) * (1.0 - psiU*c3 - 2.0*c2);
c2ddotU = 1.0/(4.0*psiU^2) * ((8.0-psiU)*c2 + 5.0*psiU*c3 - 4.0);

psiL = (lower+upper)/2*1.2;
[c2,c3] = findc2c3( psiL );
c2dotL = 1.0/(2*psiL) * (1.0 - psiL*c3 - 2.0*c2);
c2ddotL = 1.0/(4.0*psiL^2) * ((8.0-psiL)*c2 + 5.0*psiL*c3 - 4.0);

% this works, but the center point is a little left of the actual minimum,
% good starting point??
fprintf(1,'%12.5f %12.5f %12.7f %12.7f \n', psiL, psiU, c2dotL, c2dotU);
for dd = 1: 10
    psi = (psiU+psiL)*0.5;    
    [c2,c3] = findc2c3( psi );
    c2dot = 1.0/(2*psi) * (1.0 - psi*c3 - 2.0*c2);
    c2ddot = 1.0/(4.0*psi^2) * ((8.0-psi)*c2 + 5.0*psi*c3 - 4.0);
    
    % next guess determined so signs alternate on L and U
    if c2dot < 0.0 
        psiL = psi;
    end;
    if c2dot > 0.0 
        psiU = psi;
    end;
        
    fprintf(1,'%12.5f %12.5f %12.5f %14.9f %14.9f \n', psiL, psiU, c2, c2dot, c2ddot);
end
fprintf(1,'Now with Newton iteration on c2 \n');
psiU = (lower+upper)/2*0.8;
psiL = (lower+upper)/2*1.2;
psi = (psiU+psiL)*0.5;    
for dd = 1: 7
    [c2,c3] = findc2c3( psi );
    c2dot = 1.0/(2*psi) * (1.0 - psi*c3 - 2.0*c2);
    c2ddot = 1.0/(4.0*psi^2) * ((8.0-psi)*c2 + 5.0*psi*c3 - 4.0);
    
    psi = psi - c2dot/c2ddot;
     
    fprintf(1,'%12.5f %12.5f %14.9f %14.9f \n', psi, c2, c2dot, c2ddot);
end


    % ---------------------------------------------------------------------
    fprintf(1,'now for the speed tests random \n');
    pause;
    
    % test speed
    directory = 'd:\codes\library\matlab\';
    outfile = fopen(strcat(directory,'tlambert.out'), 'wt');
    ktrok = 0;
    for i=1:10000
        % if shuffle, then each will differ
        % if default, need to specify however many iterations you want and use it
        rng('default');
        a = -44000.0 + (88000.0).*rand(i,1);  % r = a + (b-a).*rand(100,1);
        b = -44000.0 + (88000.0).*rand(i,1);
        c = -44000.0 + (88000.0).*rand(i,1);
        %den = 1.0 / sqrt(a^2 + b^2 + c^2); % don't normalize
        r1 =  [a(i); b(i); c(i)];

        a = -44000.0 + (88000.0).*rand(i,1);
        b = -44000.0 + (88000.0).*rand(i,1);
        c = -44000.0 + (88000.0).*rand(i,1);     
        %den = 1.0 / sqrt(a^2 + b^2 + c^2);
        r2 =  [a(i); b(i); c(i)];

        dm = 's';
        df = 'd';
        nrev = 0;
        dtsec = 0.0 + (85000.0 - 0.0).*rand(i,1);
        %fprintf(1,' r1 %16.8f%16.8f%16.8f\n',r1 );
        %fprintf(1,' r2 %16.8f%16.8f%16.8f %11.7f \n',r2,dtsec(i) );
        [v1t,v2t,errorl] = lambertu ( r1, v1, r2, dm, df, nrev, dtwait, dtsec(i), tbi, outfile ); 
        
        if (strcmp(errorl, '      ok') ~= 0)
            ktrok = ktrok + 1;
        end     

    end
    ktrok
    
    % test cdma paper
    % case 1
    ecc = 0.2;
    a = 20000.0;
    p = a*(1.0-ecc^2);
    argp = 0.0/rad;
    [r1,v1] = coe2rv ( p, ecc, 5.0/rad, 80.0/rad, argp,     0.0, 0.0, 0.0, 0.0);
    ecc = 0.3;
    a = 30000.0;
    p = a*(1.0-ecc^2);
    argp = 30.0/rad;
    [r2,v2] = coe2rv ( p, ecc, 5.0/rad, 80.0/rad, argp,     0.0, 0.0, 0.0, 0.0);

    % case 2
    ecc = 0.6;
    a = 18000.0;
    p = a*(1.0-ecc^2);
    argp = 0.0/rad;
    [r1,v1] = coe2rv ( p, ecc, 5.0/rad, 80.0/rad, argp,     0.0, 0.0, 0.0, 0.0);
    ecc = 0.8;
    a = 35000.0;
    p = a*(1.0-ecc^2);
    argp = 120.0/rad;
    [r2,v2] = coe2rv ( p, ecc, 5.0/rad, 80.0/rad, argp,     0.0, 0.0, 0.0, 0.0);

    % case 3 no coaxial
    ecc = 0.2;
    a = 10000.0;
    p = a*(1.0-ecc^2);
    argp = 0.0/rad;
    [r1,v1] = coe2rv ( p, ecc, 5.0/rad, 80.0/rad, argp,     0.0, 0.0, 0.0, 0.0);
    ecc = 0.2;
    a = 10000.0;
    p = a*(1.0-ecc^2);
    argp = 60.0/rad;
    [r2,v2] = coe2rv ( p, ecc, 5.0/rad, 80.0/rad, argp,     0.0, 0.0, 0.0, 0.0);
    
    dm = 's';
    df = 'd';
    nrev = 0;
    dtwait = 0.0;
    
    [deltava,deltavb,dttu ] = hohmann (mag(r1)/re,mag(r2)/re,ecc,ecc,0.0,pi);      
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, dm, df, nrev, dtwait, dtsec, tbi, fid );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f  %11.6f \n',v1t, mag(v1t) );
    fprintf(1,' v2dv %16.8f%16.8f%16.8f  %11.6f \n',v2t, mag(v2t) );
    [v1dvb, v2t, errorl] = lambertb ( r1, v1, r2, dm, nrev, dtsec );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f  %11.6f \n',v1dvb, mag(v1dvb) );
    fprintf(1,' v2dv %16.8f%16.8f%16.8f  %11.6f \n',v2t, mag(v2t) );
    
    
    
%    fprintf(1,'iter       y         dtnew          psiold      psinew   psinew-psiold   dtdpsi      dtdpsi2    lower    upper      \n');
    [tbiu, tbilu] = lambgettbiu(r1, r2, 3);

    for i = 1:10
        dtsec = (200*i);

%        fprintf(1,' r1 %16.8f%16.8f%16.8f\n',r1 );
%        fprintf(1,' r2 %16.8f%16.8f%16.8f\n',r2 );
        [v1t,v2t,errorl] = lambertu ( r1, v1, r2, dm, df, nrev, dtwait, dtsec, tbiu, 1 );
%        fprintf(1,' v1t %16.8f%16.8f%16.8f\n',v1t );
%        fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    
        fprintf(1,' v1dv %16.8f%16.8f%16.8f  %11.6f\n',v1'-v1t, mag(v1'-v1));
        fprintf(1,' v2dv %16.8f%16.8f%16.8f  %11.6f\n',v2t-v2', mag(v2t-v2'));
    end
    
    
    fprintf(1,'\n-------- lambert test mark \n' );
    rad = 180.0/pi;   
    mu = 3.986004418e5;
    tusec = 806.8111238242922;
    fid = 1;
    [r1, v1] = coe2rv ( 42250.0, 0.001, 1.0/rad, 0.0/rad, 0.0, 90.0/rad,0.0,0.0/rad,0.0 );
    r1 = r1';  % needed so not to get 3 vectors in response
    [r2, v2] = coe2rv ( 40000.0, 0.001, 1.0/rad, 0.0/rad, 0.0, 87.0/rad,0.0,0.0/rad,0.0 );
    r2 = r2';  % needed so not to get 3 vectors in response
    magr1 = mag(r1);
    magr2 = mag(r2);
    % this value stays constant in all calcs, vara changes with df
    cosdeltanu = dot(r1,r2) / (magr1*magr2);  
    
    fprintf(1,'r1 %16.8f%16.8f%16.8f v1 %16.8f%16.8f%16.8f\n',r1, v1 );
    fprintf(1,'r2 %16.8f%16.8f%16.8f v2 %16.8f%16.8f%16.8f\n',r2, v2 );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1, v1, re, mu);
    fprintf(1,'ans coes p %11.4f a %11.4f e %13.9f i %13.7f raan %11.5f argp %11.5f nu %11.3f \n',...
            p,a,ecc,incl*rad,omega*rad,argp*rad, nu*rad );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2, v2, re, mu);
    fprintf(1,'ans coes p %11.4f a %11.4f e %13.9f i %13.7f raan %11.5f argp %11.5f nu %11.3f \n',...
            p,a,ecc,incl*rad,omega*rad,argp*rad, nu*rad );
    
    dtsec = 6.0*3600.0;  % sec
    dm = 's';  % 's'
    df = 'd';  % 'd'
    nrev = 0;
    dtwait = 0.0;
    [ tbi, tbil] = lambgettbiu(r1, r2, 5);
    [tbik, tbilk] = lambgettbik(r1, r2, 5);
    
    fprintf(1,'\n-------- lambertu test \n' );
    [v1t, v2t, errorl] = lambertu ( r1, v1, r2, dm, df, nrev, dtwait, dtsec, tbi, fid );
    fprintf(1,' v1t %16.8f%16.8f%16.8f\n',v1t );
    fprintf(1,' v2t %16.8f%16.8f%16.8f\n',v2t );
    [v1tb, v2tb, errorl] = lambertb ( r1, v1, r2, dm, nrev, dtsec );
    fprintf(1,' v1tb %16.8f%16.8f%16.8f\n',v1tb );
    fprintf(1,' v2tb %16.8f%16.8f%16.8f\n',v2tb );
    [v1tk, v2tk, errorl] = lambertk ( r1, v1, r2, dm, df, nrev, dtwait, dtsec, tbik, fid );
    fprintf(1,'\n v1tk %16.8f%16.8f%16.8f\n',v1tk );
    fprintf(1,' v2tk %16.8f%16.8f%16.8f\n\n',v2tk );

    dv1 = v1'-v1t;
    mag(dv1);
    dv2 = v2'-v2t;
    mag(dv2);
    fprintf(1,' dv1 %16.8f%16.8f%16.8f\n',dv1 );
    fprintf(1,' dv2 %16.8f%16.8f%16.8f\n',dv2 );
    

     fprintf(1,'-------------------- problem ex 6-3 \n');
     rinit  = (42250.0)/re;
     rfinal = (40000.0)/re;
     einit  = 0.0;
     efinal = 0.0;
     nuinit = 0.0/rad;
     nutran = 160.0/rad;

     fprintf(1,'initial position \n');
     fprintf(1,' rinit  %11.7f %11.7f km \n',rinit, rinit*re);
     fprintf(1,' rfinal %11.7f %11.7f km \n',rfinal, rfinal*re);
     fprintf(1,' einit   %11.7f \n',einit);
     fprintf(1,' efinal  %11.7f \n',efinal);
     fprintf(1,' nuinit  %11.7f deg \n',nuinit * rad);
     fprintf(1,' nutran %11.7f deg \n',nutran * rad);

     [deltava,deltavb,dttu,etran,atran, vtrana, vtranb ] = onetang(rinit,rfinal,einit,efinal,nuinit,nutran);

     constastro;
     fprintf(1,'one tangent answers \n');
     fprintf(1,' deltava  %11.7f  %11.7f km/s \n',deltava, deltava*velkmps );
     fprintf(1,' deltavb  %11.7f  %11.7f km/s \n',deltavb, deltavb*velkmps );
     fprintf(1,' deltav  %11.7f %11.7f   km/s \n',deltavb + deltava, (deltava + deltavb)*velkmps );
     fprintf(1,' dttu  %11.7f tu %11.7f min \n',dttu,dttu*tumin);

    
    