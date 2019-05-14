%
% convert orbital elments and find nu
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
% jdut1=jday(2014, 7, 26, 0, 0, 0.0);
% nu = lon2nu (jdut1, 7.020438698/rad, 0.070273056/rad, 19.90450011/rad, 352.5056022/rad);
%
function nu = lon2nu( jdut1, lon, incl, raan, argp)
     rad = 180/pi;
     twopi = 2.0 * pi;

 %    fprintf(' jd %16.8f lon %11.5f  incl %11.5f raan %11.5f argp %11.5f \n',jdut1, lon*rad, incl*rad, raan*rad, argp*rad );
 %    need to use their GMST calculation
     ed = jdut1 +0.0 - 2451544.5;  % elapsed days from 1 jan 2000 0 hrs
     gmst = 99.96779469 + 360.9856473662860 * ed + 0.29079e-12 * ed * ed;  % deg
     deg2rad    = pi/180.0;
     gmst = rem( gmst*deg2rad,2.0*pi );

     % ------------------------ check quadrants --------------------
     if ( gmst < 0.0 )
         gmst = gmst + 2.0*pi;
     end

     lambdau = gmst + lon - raan;

     % make sure lambdau is 0 to 360 deg
     if lambdau < 0.0  
         lambdau = lambdau + 2.0*pi;
     end
     if lambdau > twopi
         lambdau = lambdau - 2.0*pi;
     end

     arglat = atan(tan(lambdau) / cos(incl));  
     % find nu     
     if (lambdau >= 0.5*pi) && (lambdau <  1.5*pi) 
         arglat = arglat + pi; 
     end
      
     temp = arglat - argp;
     
     nu = temp;
% fprintf(' %11.5f %11.5f %11.5f %11.5f  %16.10f ',lambdau*rad, argp*rad, lon*rad, gmst*rad, temp*rad);
 fprintf(' lu %11.5f argp %11.5f lon %11.5f gmst %11.5f arglat %11.5f nu %16.10f ',lambdau*rad, argp*rad, lon*rad, gmst*rad, arglat*rad, temp*rad);
    
 %     fprintf(' nu = %11.5f deg \n',nu*rad);
end