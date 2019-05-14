%     -----------------------------------------------------------------
%
%                              ex6_7.m
%
%  this file demonstrates example 6-7.
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
      vkmps = 7.905365719014;

      fprintf(1,'-------------------- problem ex 6-7 \n');
      rinit  = (re + 191.0)/re;
      rfinal = (re + 35780.0)/re;
      fprintf(1,'from radius %11.5f to %11.5f ---------------------- \n', rinit*re, rfinal*re);
      einit  = 0.0;
      efinal = 0.0;
      iinit = 28.5/rad;
      ifinal=  0.0/rad;
      fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
      deltai = -ifinal + iinit;
      nuinit =  0.0/rad;
      nufinal=  180.0/rad;
         
      [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f %11.7f km/s \n',deltava, deltava*vkmps );
      fprintf(1,' deltavb  %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps );
      fprintf(1,' deltai1  %11.7f  %11.7f  \n',deltai1 * rad, deltai2 * rad);
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

      atran = (rinit+rfinal)*re/2.0;
      etran = re*(rfinal-rinit)/(2*atran);
      ptran = atran*(1.0-etran^2);
      %                   a      e       i               node         argp    nu      ci arglat ce truelon ee lonper             
      [r1,v1] = coe2rv(rinit*re,0.0, 28.5/rad,         45.0/rad,   0.0/rad,   0.0/rad,  0.0/rad, 0.0,0.0);
      [r2,v2] = coe2rv(ptran, etran, 28.5/rad-deltai1, 45.0/rad,   0.0/rad,   0.0/rad,  0.0,    0.0,    0.0);
      [r3,v3] = coe2rv(ptran, etran, 28.5/rad-deltai1, 45.0/rad,   0.0/rad, 180.0/rad,  0.0,    0.0,    0.0);
      [r4,v4] = coe2rv(rfinal*re, 0.0, ifinal,          0.0/rad,   0.0/rad,   0.0/rad,  0.0,  225.0/rad, 0.0);
      
      fprintf(1,' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v4-v3, mag(v4-v3) );


      % ===================================================================
      rinit  = (re + 222.0)/re;
      rfinal = (re + 35786.0)/re;
      einit  = 0.0;
      efinal = 0.0;
      iinit = 28.5/rad;
      ifinal= 89.0/rad;
      deltai = ifinal - iinit;
      nuinit =  0.0/rad;
      nufinal=  180.0/rad;
      fprintf(1,'from radius %11.5f to %11.5f GEO Polar ---------------------- \n', rinit*re, rfinal*re);
      fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
         
      [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f %11.7f km/s \n',deltava, deltava*vkmps );
      fprintf(1,' deltavb  %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps );
      fprintf(1,' deltai1  %11.7f  %11.7f  \n',deltai1 * rad, deltai2 * rad);
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

      atran = (rinit+rfinal)*re/2.0;
      etran = re*(rfinal-rinit)/(2*atran);
      ptran = atran*(1.0-etran^2);
      %                   a      e       i               node         argp    nu      ci arglat ce truelon ee lonper             
      [r1,v1] = coe2rv(rinit*re,0.0, 28.5/rad,         45.0/rad,   0.0/rad,   0.0/rad,    0.0/rad, 0.0,0.0);
      [r2,v2] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad,   0.0/rad,    0.0,    0.0,    0.0);
      [r3,v3] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad, 180.0/rad,    0.0,    0.0,    0.0);
      [r4,v4] = coe2rv(rfinal*re, 0.0, ifinal,         45.0/rad,   0.0/rad,   0.0/rad,  180.0,    0.0, 0.0);
      
      fprintf(1,' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v4-v3, mag(v4-v3) );
      
      % ===================================================================
      rinit  = (re + 222.0)/re;
      rfinal = (re + 35786.0)/re;
      einit  = 0.0;
      efinal = 0.0;
      iinit = 28.5/rad;
      ifinal= 0.0/rad;
      deltai = ifinal - iinit;
      nuinit =  0.0/rad;
      nufinal=  180.0/rad;
      fprintf(1,'from radius %11.5f to %11.5f GEO ---------------------- \n', rinit*re, rfinal*re);
      fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
         
      [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f %11.7f km/s \n',deltava, deltava*vkmps );
      fprintf(1,' deltavb  %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps );
      fprintf(1,' deltai1  %11.7f  %11.7f  \n',deltai1 * rad, deltai2 * rad);
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

      atran = (rinit+rfinal)*re/2.0;
      etran = re*(rfinal-rinit)/(2*atran);
      ptran = atran*(1.0-etran^2);
      %                   a      e       i               node         argp    nu       arglat  truelon lonper             
      [r1,v1] = coe2rv(rinit*re,0.0, 28.5/rad,         45.0/rad,   0.0/rad, 0.0/rad,   0.0/rad, 0.0,0.0);
      [r2,v2] = coe2rv(ptran, etran, 28.5/rad-deltai1, 45.0/rad,   0.0/rad, 0.0/rad,    0.0,    0.0,    0.0);
      [r3,v3] = coe2rv(ptran, etran, 28.5/rad-deltai1, 45.0/rad,   0.0/rad, 180.0/rad,  0.0,    0.0,    0.0);
      [r4,v4] = coe2rv(rfinal*re, 0.0, ifinal,          0.0/rad,   0.0/rad,   0.0/rad,  0.0,  225.0/rad, 0.0);
      
      fprintf(1,' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v4-v3, mag(v4-v3) );
      
      % ===================================================================
      rinit  = (re + 222.0)/re;
      rfinal = (re + 15392.0)/re;
      einit  = 0.0;
      efinal = 0.64887;
      iinit = 28.5/rad;
      ifinal= 63.4/rad;
      deltai = ifinal - iinit;
      nuinit =  0.0/rad;
      nufinal=  180.0/rad;
      fprintf(1,'from radius %11.5f to %11.5f Mod HEO ---------------------- \n', rinit*re, rfinal*re);
      fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
         
      [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f %11.7f km/s \n',deltava, deltava*vkmps );
      fprintf(1,' deltavb  %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps );
      fprintf(1,' deltai1  %11.7f  %11.7f  \n',deltai1 * rad, deltai2 * rad);
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

      atran = (rinit+rfinal)*re/2.0;
      etran = re*(rfinal-rinit)/(2*atran);
      ptran = atran*(1.0-etran^2);
      %                   a      e       i               node         argp    nu      ci arglat ce truelon ee lonper             
      [r1,v1] = coe2rv(rinit*re,0.0, 28.5/rad,         45.0/rad,   0.0/rad,   0.0/rad,    0.0/rad, 0.0,0.0);
      [r2,v2] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad,   0.0/rad,    0.0,    0.0,    0.0);
      [r3,v3] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad, 180.0/rad,    0.0,    0.0,    0.0);
      [r4,v4] = coe2rv(rfinal*re, efinal, ifinal,      45.0/rad, 270.0/rad, 180.0/rad,    0.0,    0.0, 0.0);
      
      fprintf(1,' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v4-v3, mag(v4-v3) );
      
      % ===================================================================
      rinit  = (re + 222.0)/re;
      rfinal = (re + 35786.0)/re;
      einit  = 0.0;
      efinal = 0.0;
      iinit = 28.5/rad;
      ifinal= 49.3/rad;
      deltai = ifinal - iinit;
      nuinit =  0.0/rad;
      nufinal=  180.0/rad;
      fprintf(1,'from radius %11.5f to %11.5f Quasi-Zen ---------------------- \n', rinit*re, rfinal*re);
      fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
         
      [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f %11.7f km/s \n',deltava, deltava*vkmps );
      fprintf(1,' deltavb  %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps );
      fprintf(1,' deltai1  %11.7f  %11.7f  \n',deltai1 * rad, deltai2 * rad);
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

      atran = (rinit+rfinal)*re/2.0;
      etran = re*(rfinal-rinit)/(2*atran);
      ptran = atran*(1.0-etran^2);
      %                   a      e       i               node         argp    nu       arglat  truelon lonper             
      [r1,v1] = coe2rv(rinit*re,0.0, 28.5/rad,         45.0/rad,   0.0/rad,   0.0/rad,    0.0/rad, 0.0,0.0);
      [r2,v2] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad,   0.0/rad,    0.0,    0.0,    0.0);
      [r3,v3] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad, 180.0/rad,    0.0,    0.0,    0.0);
      [r4,v4] = coe2rv(rfinal*re, 0.0, ifinal,         45.0/rad,   0.0/rad,   0.0/rad,  180.0,    0.0, 0.0);
      
      fprintf(1,' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v4-v3, mag(v4-v3) );

      % ===================================================================
      rinit  = (re + 222.0)/re;
      rfinal = (re + 37028.0)/re;
      einit  = 0.0;
      efinal = 0.0;
      iinit = 28.5/rad;
      ifinal= 0.0/rad;
      deltai = ifinal - iinit;
      nuinit =  0.0/rad;
      nufinal=  180.0/rad;
      fprintf(1,'from radius %11.5f to %11.5f Super GEO ---------------------- \n', rinit*re, rfinal*re);
      fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
         
      [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f %11.7f km/s \n',deltava, deltava*vkmps );
      fprintf(1,' deltavb  %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps );
      fprintf(1,' deltai1  %11.7f  %11.7f  \n',deltai1 * rad, deltai2 * rad);
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

      atran = (rinit+rfinal)*re/2.0;
      etran = re*(rfinal-rinit)/(2*atran);
      ptran = atran*(1.0-etran^2);
      %                   a      e       i               node         argp    nu       arglat  truelon lonper             
      [r1,v1] = coe2rv(rinit*re,0.0, 28.5/rad,         45.0/rad,   0.0/rad,   0.0/rad,    0.0/rad, 0.0,0.0);
      [r2,v2] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad,   0.0/rad,    0.0,    0.0,    0.0);
      [r3,v3] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad, 180.0/rad,    0.0,    0.0,    0.0);
      [r4,v4] = coe2rv(rfinal*re, 0.0, ifinal,          0.0/rad,   0.0/rad,   0.0/rad,  0.0,  225.0/rad, 0.0);
      
      fprintf(1,' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v4-v3, mag(v4-v3) );
    
      % ===================================================================
      rinit  = (re + 222.0)/re;
      rfinal = (re + 35786.0)/re;
      einit  = 0.0;
      efinal = 0.268;
      iinit = 28.5/rad;
      ifinal= 63.4/rad;
      deltai = ifinal - iinit;
      nuinit =  0.0/rad;
      nufinal=  180.0/rad;
      fprintf(1,'from radius %11.5f to %11.5f Sirius ---------------------- \n', rinit*re, rfinal*re);
      fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
         
      [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f %11.7f km/s \n',deltava, deltava*vkmps );
      fprintf(1,' deltavb  %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps );
      fprintf(1,' deltai1  %11.7f  %11.7f  \n',deltai1 * rad, deltai2 * rad);
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

      atran = (rinit+rfinal)*re/2.0;
      etran = re*(rfinal-rinit)/(2*atran);
      ptran = atran*(1.0-etran^2);
      %                   a      e       i               node         argp    nu       arglat  truelon lonper             
      [r1,v1] = coe2rv(rinit*re,0.0, 28.5/rad,         45.0/rad,   0.0/rad,   0.0/rad,    0.0/rad, 0.0,0.0);
      [r2,v2] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad,   0.0/rad,    0.0,    0.0,    0.0);
      [r3,v3] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad, 180.0/rad,    0.0,    0.0,    0.0);
      [r4,v4] = coe2rv(rfinal*re, efinal, ifinal,      45.0/rad, 270.0/rad, 180.0/rad,    0.0,    0.0, 0.0);
      
      fprintf(1,' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v4-v3, mag(v4-v3) );
      
      % ===================================================================
      rinit  = (re + 222.0)/re;
      rfinal = (re + 35786.0)/re;
      einit  = 0.0;
      efinal = 0.700381;
      iinit = 28.5/rad;
      ifinal= 63.4/rad;
      deltai = ifinal - iinit;
      nuinit =  0.0/rad;
      nufinal=  180.0/rad;
      fprintf(1,'from radius %11.5f to %11.5f Tundra-like ---------------------- \n', rinit*re, rfinal*re);
      fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
         
      [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f %11.7f km/s \n',deltava, deltava*vkmps );
      fprintf(1,' deltavb  %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps );
      fprintf(1,' deltai1  %11.7f  %11.7f  \n',deltai1 * rad, deltai2 * rad);
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

      atran = (rinit+rfinal)*re/2.0;
      etran = re*(rfinal-rinit)/(2*atran);
      ptran = atran*(1.0-etran^2);
      %                   a      e       i               node         argp    nu       arglat  truelon lonper             
      [r1,v1] = coe2rv(rinit*re,0.0, 28.5/rad,         45.0/rad,   0.0/rad,   0.0/rad,    0.0/rad, 0.0,0.0);
      [r2,v2] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad,   0.0/rad,    0.0,    0.0,    0.0);
      [r3,v3] = coe2rv(ptran, etran, 28.5/rad+deltai1, 45.0/rad,   0.0/rad, 180.0/rad,    0.0,    0.0,    0.0);
      [r4,v4] = coe2rv(rfinal*re, efinal, ifinal,      45.0/rad, 270.0/rad, 180.0/rad,    0.0,    0.0, 0.0);
      
      fprintf(1,' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v2-v1, mag(v2-v1) );
      fprintf(1,' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n',v4-v3, mag(v4-v3) );
      