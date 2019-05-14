%     -----------------------------------------------------------------
%
%                              Ex11_4.m
%
%  this file demonstrates example 11-4.
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
%            30 aug 11  david vallado
%                         original
%  changes :
%            22 aug 11 david vallado
%                         original baseline
%
%     *****************************************************************

        constmath;
        constastro;
        j2 = 0.00108263;
        
        % --------  repeat gt calculations
        a = 6570.3358; % km
        ecc = 0.006301;
        
        incl = 45.00/rad;

        p = a * (1.0 - ecc*ecc);
        nanom = sqrt( mu/(a*a*a) );
   fprintf(1,'p %11.7f  %11.7f \n',p, p/re);
   
        n = nanom;
   fprintf(1,'n %11.7f  %11.7f \n',n, n/rad);
        
        raanrate = -1.5*j2*re^2*n*cos(incl)/(p*p);
   fprintf(1,'raanrate %11.7f  %11.7f \n',raanrate, raanrate*180*86400/pi);  % raanrate/omegaearthradptu

        deltatoa = (3.0*pi/(n*a)) * (1.0 + 0.5*j2*(re/a)^2*(4.0*(cos(incl))^2 - 1.0));
   fprintf(1,'deltatoa %11.7f  %11.7f \n',deltatoa, deltatoa*re/tusec);
        
        deltaOdotoa = -3.5*raanrate/a;
   fprintf(1,'deltaOdotoa %d  %11.11f \n',deltaOdotoa, deltaOdotoa*re/tusec);
        
        deltapoi = 12.0*pi/n * j2 * (re/a)^2 * sin(2.0*incl);
   fprintf(1,'deltapoi %11.7f  %11.7f \n',deltapoi, deltapoi/tusec);
        
        deltaraanoi = -raanrate * tan(incl);
   fprintf(1,'deltaraanoi %11.7f  %11.7f \n',deltaraanoi, deltaraanoi/tusec);

        pnodal = 2.0 * pi / n;
   fprintf(1,'pnodal %11.7f s %11.7f min \n',pnodal, pnodal/60.0);
        anomper = pnodal * (1.0 / (1.0 + 0.75 * j2 * (re/p)^2 * (sqrt(1.0 - ecc*ecc) * (2.0 - 3.0 * sin(incl)^2) + (4.0 - 5.0 * sin(incl)^2 ) ) ) );
   fprintf(1,'anomper %11.7f s %11.7f min \n',anomper, anomper/60.0);

        dellon = (omegaearth - raanrate)*anomper;
   fprintf(1,'dellon %11.7f  %11.7f \n',dellon, dellon*re);
        dlpa = re * (omegaearth - raanrate) * deltatoa - deltaOdotoa * re * anomper;
   fprintf(1,'dlpa %11.7f  %11.7f \n',dlpa, dlpa);
        dlpi = re * (omegaearth - raanrate) * deltapoi - deltaraanoi * re * anomper;
   fprintf(1,'dlpi %11.7f  %11.7f \n',dlpi, dlpi);

% osculating values from stk
%2 Apr 2000 12:00:00.000  6570.3400000 0.00630100  45.00000  321.21600  69.629    0.000   5300.206   
%2 Apr 2000 13:28:20.000  6565.1860071 0.00579447  44.99776  320.82382  70.494  359.550   5293.971  
% mean values from stk


        dadt = -0.174562 / anomper;
        didt = 0.000132 * pi / (180 * anomper);
%
%        dadt = -1.59082 / anomper;
%        didt =  0.01188 * pi / (180 * anomper);

   fprintf(1,'dadt %11.7f km/rev %d km/s \n',dadt*anomper, dadt);
   fprintf(1,'didt %11.7f deg/rev %d deg/sec \n',didt*anomper*180/pi, didt);
        
        k2 = 1.0/anomper * (dlpa * dadt + dlpi * didt);
        k1 = sqrt(2.0 * k2 * (-50 - 50) );
   fprintf(1,'k2 %d  %11.7f \n',k2, k2*tusec/re);
   fprintf(1,'k1 %d  %11.7f \n',k1, k1*tusec/re);
   
        da = k1 * anomper / dlpa;
   fprintf(1,'da %11.7f  %11.7f \n',da, da);
        
        tdrift = -k1/k2;
   fprintf(1,'tdrift %11.7f sec  %11.7f min %11.7f day \n',tdrift, tdrift/60.0, tdrift/86400.0);
        deltav = n * 0.5 * da;
   fprintf(1,'deltav %11.7f km/s %11.7f m/s \n',deltav, deltav*1000);
       

      fprintf(1,'change frequency deltav %11.7f m/s timing  %11.7f min \n',2*deltav*1000, 2*tdrift/60);

      x = x/0;
        
        

