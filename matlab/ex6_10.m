%     -----------------------------------------------------------------
%
%                              ex6_10.m
%
%  this file demonstrates example 6-10.
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

      constastro;
      
      rad = 180.0 / pi;
      re = 6378.137;  
      fprintf(1,'-------------------- problem ex 6-10 \n');
      aint = 7143.51 / re;
      atgt = 42159.4855 / re; 
      iint = 28.5 / rad;
      itgt = 0.0 / rad;
      deltai = itgt - iint;
      nodeint = 45.0 / rad;
      arglatint = 15.0 / rad;
      phasenew = pi - arglatint;
      truelon = 200.0 / rad;
      ktgt = 0;
      kint = 1;

      fprintf(1,'combined maneuver \n');
      
      [ttrans,tphase,dvphase,dvtrans1,dvtrans2,aphase ] = noncoplr(phasenew,aint,atgt,ktgt,kint,arglatint,nodeint,truelon,deltai);

      fprintf(1,' ttrans  %11.7f \n',ttrans );
      fprintf(1,' tphase  %11.7f \n',tphase );
      fprintf(1,' dvphase  %11.7f  %11.7f \n',dvphase, dvphase*velkmps );
      fprintf(1,' dvtrans1  %11.7f  %11.7f \n',dvtrans1, dvtrans1*velkmps );
      fprintf(1,' dvtrans2  %11.7f  %11.7f \n',dvtrans2, dvtrans2*velkmps );
      fprintf(1,' aphase  %11.7f \n',aphase );