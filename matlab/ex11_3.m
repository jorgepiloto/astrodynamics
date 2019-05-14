%     -----------------------------------------------------------------
%
%                              Ex11_3.m
%
%  this file demonstrates example 11-3.
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
%            30 dec 12  david vallado
%                         original
%  changes :
%            30 dec 12 david vallado
%                         original baseline
%
%     *****************************************************************

        constmath;
        constastro;
        j2 = 0.00108263;
        
        % --------  repeat gt calculations
        alt = 160.0; % km
        ecc = 0.006301;
        incl = 45.0/rad;
        
        latgd = 41.52/rad;
        lon = 12.37/rad;
        salt = 0.152; %km
        
        etafov = 20.0/rad;
        etactr = 15.0/rad;
        
        sinlat      = sin( latgd );
        % ------  find rdel and rk components of site vector  ---------
        cearth= re / sqrt( 1.0 - ( eccearthsqrd*sinlat*sinlat ) );
        rdel  = ( cearth + salt )*cos( latgd );
        rk    = ( (1.0-eccearthsqrd)*cearth + salt )*sinlat;

        % ---------------  find site position vector  -----------------
        rs(1) = rdel * cos( lon );
        rs(2) = rdel * sin( lon );
        rs(3) = rk;

        r = mag(rs);
        rp = r + alt;
      fprintf(1,'rdel %11.7f  rk %11.7f r %11.7f  rpsite %11.7f km \n', rdel, rk, r, rp);
      
        r = re + alt;

        fovmin = 0.5*etafov + etactr;
        gamma   = pi - asin( r*sin(fovmin)/re );  % use larger angle
   fprintf(1,'gamma %11.7f gamma %11.7f rho %11.7f fovmin %11.7f \n', gamma*rad, (pi-gamma)*rad, r, fovmin);
        
        rho  = re*cos( gamma ) + r*cos(fovmin);
   fprintf(1,'fovmin %11.7f gamma %11.7f gamma %11.7f rho %11.7f  \n',fovmin*rad, gamma*rad, (pi-gamma)*rad, rho);

        lambda  = asin(rho * sin(fovmin)/re);
   fprintf(1,'lambda %11.7f  %11.7f km  \n',lambda*rad, lambda*re);
   
        revpday = 16.4;
        n = 16.4*2*pi/86400.0;  % rad/s
        a = (mu*(1/n)^2)^0.33333;
        ecc = (a-rp)/a;
   fprintf(1,'n %11.7f r/d %11.7f rad/s  a %11.7f km ecc %11.7f  \n',revpday, n, a, ecc);
        
   fprintf(1,'\n\n now start the deisgn process letting ecc = 0.001 \n\n');
        ecc = 0.001;
        a = rp/(1.0-ecc);
        n = sqrt(mu/(a*a*a));
        omega = n/omegaearth;
   fprintf(1,'n %11.7f rad/s  a %11.7f km ecc %11.7f omega  %11.7f \n', n, a, ecc, omega);
        
        a = 6535.4713;
        ecc = 0.001;
        pii = 16.387;
        qi = 1;
        p = a*(1.0 - ecc*ecc);
        n = sqrt(mu/(a*a*a));
        raandot = -1.5*j2*(re/p)^2*n*cos(incl);
        nn = pii/qi * (omegaearth - raandot) * (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
   
        kepperiod = (2.0 * pi) / n;
        nodalperiodG = (2.0 * pi) / (omegaearth - raandot);
        nodalperiod = kepperiod / (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
        
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/60, nodalperiod/60, kepperiod/60, nn);
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/tusec, nodalperiod/tusec, kepperiod/tusec, nn);
        

        pii = 16.0;
        qi = 1;
        p = a*(1.0 - ecc*ecc);
        n = sqrt(mu/(a*a*a));
        raandot = -1.5*j2*(re/p)^2*n*cos(incl);
        nn = pii/qi * (omegaearth - raandot) * (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
   
        kepperiod = (2.0 * pi) / n;
        nodalperiodG = (2.0 * pi) / (omegaearth - raandot);
        nodalperiod = kepperiod / (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
        
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/60, nodalperiod/60, kepperiod/60, nn);
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/tusec, nodalperiod/tusec, kepperiod/tusec, nn);

        pii = 32;
        qi = 2;
        p = a*(1.0 - ecc*ecc);
        n = sqrt(mu/(a*a*a));
        raandot = -1.5*j2*(re/p)^2*n*cos(incl);
        nn = pii/qi * (omegaearth - raandot) * (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
   
        kepperiod = (2.0 * pi) / n;
        nodalperiodG = (2.0 * pi) / (omegaearth - raandot);
        nodalperiod = kepperiod / (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
        
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/60, nodalperiod/60, kepperiod/60, nn);
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/tusec, nodalperiod/tusec, kepperiod/tusec, nn);

        pii = 49;
        qi = 3;
        p = a*(1.0 - ecc*ecc);
        n = sqrt(mu/(a*a*a));
        raandot = -1.5*j2*(re/p)^2*n*cos(incl);
        nn = pii/qi * (omegaearth - raandot) * (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
   
        kepperiod = (2.0 * pi) / n;
        nodalperiodG = (2.0 * pi) / (omegaearth - raandot);
        nodalperiod = kepperiod / (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
        
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/60, nodalperiod/60, kepperiod/60, nn);
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/tusec, nodalperiod/tusec, kepperiod/tusec, nn);

        pii = 82;
        qi = 5;
        p = a*(1.0 - ecc*ecc);
        n = sqrt(mu/(a*a*a));
        raandot = -1.5*j2*(re/p)^2*n*cos(incl);
        nn = pii/qi * (omegaearth - raandot) * (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
   
        kepperiod = (2.0 * pi) / n;
        nodalperiodG = (2.0 * pi) / (omegaearth - raandot);
        nodalperiod = kepperiod / (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
        
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/60, nodalperiod/60, kepperiod/60, nn);
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/tusec, nodalperiod/tusec, kepperiod/tusec, nn);
   
        pii = 213;
        qi = 33;
        p = a*(1.0 - ecc*ecc);
        n = sqrt(mu/(a*a*a));
        raandot = -1.5*j2*(re/p)^2*n*cos(incl);
        nn = pii/qi * (omegaearth - raandot) * (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
   
        kepperiod = (2.0 * pi) / n;
        nodalperiodG = (2.0 * pi) / (omegaearth - raandot);
        nodalperiod = kepperiod / (1 - 1.5*j2*(re/a)^2*(3.0 - 4.0*sin(incl)*sin(incl)));
        
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/60, nodalperiod/60, kepperiod/60, nn);
   fprintf(1,'p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n', pii, qi, nodalperiodG/tusec, nodalperiod/tusec, kepperiod/tusec, nn);
   
        
        
        

