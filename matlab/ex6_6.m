%     -----------------------------------------------------------------
%
%                              ex6_6.m
%
%  this file demonstrates example 6-6.
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
%            25 nov 08  david vallado
%                         original
%  changes :
%            25 nov 08  david vallado
%                         original baseline
%
%     *****************************************************************

      rad = 180.0 / pi;
      re = 6378.137;  
      mu = 1.0;  % canonical

      fprintf(1,'-------------------- problem ex 6-6 \n');
      iinit= 55.0 / rad;
      ifinal = 40.0 / rad;
      ecc = 0.0;
      deltaomega = 45.0 / rad;
      vinit = 5.892311 / 7.905365719;
      fpa = 0.0 / rad;

      [deltav] = iandnode(iinit,deltaomega,ifinal,vinit,fpa);

      fprintf(1,'inclination and node changes \n');
      fprintf(1,' deltav %11.7f %11.7f km/s \n',deltav, deltav*7.905365719 );

      
      [r1,v1] = coe2rv(11480.649,0.0,55.0/rad,45.0/rad,0.0/rad,330.0/rad,128.9041397/rad,0.0,0.0);
      [r2,v2] = coe2rv(11480.649,0.0,40.0/rad,90.0/rad,0.0/rad,330.0/rad,97.3803453/rad,0.0,0.0);
      
      fprintf(1,' form Dv vectors and look at magnitude difference \n comp %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
          