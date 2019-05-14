% 
% convert orbital elments and find lon
%  input all angles in rad
%  output is in rad
%  dav 6 may 2011
%
% jdut1=jday(2011,3,22,12,30,0);
% incl = 0.04801366433780/rad; 
% raan =  330.034263262230/rad; 
% argp =  91.3766761083524/rad; 
% lon =  -7.2935047164727/rad;
%
function lon = nu2lon( jdut1, nu, incl, raan, argp)
     rad = 180/pi;
     twopi = 2.0 * pi;

 %    fprintf(' jd %16.8f lon %11.5f  incl %11.5f raan %11.5f argp %11.5f \n',jdut1, lon*rad, incl*rad, raan*rad, argp*rad );
 %    need to use their GMST calculation
     ed= jdut1 - 2451544.5;  % elapsed days from 1 jan 2000 0 hrs
     gmst = 99.96779469 + 360.9856473662860 * ed + 0.29079e-12 * ed * ed;  % deg
     deg2rad    = pi/180.0;
     gmst = rem( gmst*deg2rad,2.0*pi );

     % ------------------------ check quadrants --------------------
     if ( gmst < 0.0 )
         gmst = gmst + 2.0*pi;
     end
     
     arglat = nu + argp;
     % make sure lambdau is 0 to 360 deg
     if arglat < 0.0  
         arglat = arglat + 2.0*pi;
     end
     if arglat > twopi
         arglat = arglat - 2.0*pi;
     end
    
     lambdau = atan(tan(arglat) * cos(incl));  
     if (arglat >= 0.5*pi) && (arglat <  1.5*pi) 
         lambdau = lambdau + pi; 
     end

     temp = lambdau - gmst + raan;
% fprintf(' xx %11.5f  %11.5f %11.5f %11.5f ', lambdau*rad, arglat*rad, raan*rad, temp*rad );
     % make sure lambdau is 0 to 360 deg
     temp = rem( temp, 2.0*pi );

     lon = temp;
% fprintf(' %11.5f  %11.5f ', nu*rad, raan*rad );
    
 %     fprintf(' nu = %11.5f deg \n',nu*rad);
end