%     -----------------------------------------------------------------
%
%                              ex6_5.m
%
%  this file demonstrates example 6-5.
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

      fprintf(1,'-------------------- problem ex 6-5 \n');
      rad = 180.0 / pi;
      re = 6378.137;  
      mu = 1.0;  % canonical

      iinit= 55.0 / rad;
      ecc = 0.0;
      deltaomega = 45.0 / rad;
      vinit = 5.892311 / 7.905365719;
      fpa = 0.0 / rad;
      incl = 55.0 / rad;

      [ifinal,deltav ] = nodeonly(iinit,ecc,deltaomega,vinit,fpa,incl);

      fprintf(1,'node only changes \n');
      fprintf(1,' ifinal %11.7f \n',ifinal*rad );
      fprintf(1,' deltav %11.7f %11.7f \n',deltav, deltav*7.905365719 );

      %                                           node       argp      nu            ci           ce   ee
      [r1,v1] = coe2rv(11480.649, 0.0, 55.0/rad, 45.0/rad, 30.0/rad, 330.0/rad, 103.3647275/rad, 0.0, 0.0);
      [r2,v2] = coe2rv(11480.649, 0.0, 55.0/rad, 90.0/rad, 30.0/rad, 330.0/rad,  76.6352725/rad, 0.0, 0.0);
      
      fprintf(1,' comp %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );


      iinit= 90.0 / rad;
      ecc = 0.0;
      deltaomega = 0.98564736 / rad;  % rad/sec
      vinit = sqrt(398600.4418/42164)/7.905365719;  % er/tu
      fpa = 0.0 / rad;
      incl = 90.0 / rad;

      [ifinal,deltav ] = nodeonly(iinit,ecc,deltaomega,vinit,fpa,incl);

      fprintf(1,'node only changes \n');
      fprintf(1,' ifinal %11.7f \n',ifinal*rad );
      fprintf(1,' deltav %11.7f %11.7f km/s \n',deltav, deltav*7.905365719);
      
      %                                       node                argp      nu          ci    ce   ee
      [r1,v1] = coe2rv(42164,0.0, 90.0/rad, 0.0/rad,            0.0/rad, 0.0/rad, 90.0/rad, 0.0, 0.0);
      [r2,v2] = coe2rv(42164,0.0, 90.0/rad, 0.0/rad+deltaomega, 0.0/rad, 0.0/rad, 90.0/rad, 0.0, 0.0);
      
      fprintf(1,' manv %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );


      