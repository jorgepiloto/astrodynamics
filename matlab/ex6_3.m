%     -----------------------------------------------------------------
%
%                              ex6_3.m
%
%  this file demonstrates example 6-3.
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
      fprintf(1,'-------------------- problem ex 6-3 \n');
      rinit  = (re + 191.3411)/re;
      rfinal = (re + 35781.34857)/re;
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
       
      % ellip equatorial
      p = re* atran*(1.0-etran*etran);
      [r1,v1] = coe2rv(p,  etran,  0.0/rad, 0.0/rad, 0.0/rad, nutran,     0.0/rad, 0.0,  0.0/rad);
      % p = a for circular, cir equatorial 
      p = re + 35781.34857;
      [r2,v2] = coe2rv(p , efinal, 0.0/rad, 0.0/rad, 0.0/rad, 0.0/rad,     0.0/rad,  nutran,  0.0);
      fprintf(1,' comp %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' r1  %11.7f %11.7f %11.7f  km v %11.7f \n',r1, mag(v1) );
      fprintf(1,' r2  %11.7f %11.7f %11.7f  km v %11.7f \n',r2, mag(v2) );
      
      