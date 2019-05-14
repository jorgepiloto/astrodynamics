%     -----------------------------------------------------------------
%
%                              Ex3_1415a.m
%
%  this file demonstrates example 3-1415a - just the teme to ecef conversion
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (h)               email dvallado@msn.com
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            10 dec 07  david vallado
%                         original
%  changes :
%            10 dec 07  david vallado
%                         original baseline
%
%     *****************************************************************

%           constastro;

%    book example
           recef = [-1033.4793830;  7901.2952754;  6380.3565958 ];
           vecef = [-3.225636520;  -2.872451450;   5.531924446 ];
           aecef = [0.001;0.002;0.003];

           year=2004;
           mon = 4;
           day = 6;
           hr =  7;
           min= 51;
           sec= 28.386009;

           dut1 = -0.4399619;
           dat  = 32;
           xp   = -0.140682;  % "
           yp   =  0.333309;
           lod  =  0.0015563;
           ddpsi = -0.052195;  % "
           ddeps = -0.003875;
           timezone = 0;
           order = 106;
           eqeterms = 2; % use the extra ee terms in j2000

           [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
           = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
           fprintf(1,'ut1 %8.6f tut1 %16.12f jdut1 %18.11f ',ut1,tut1,jdut1+jdut1frac );
           [h,m,s] = sec2hms( ut1 );
           fprintf(1,'hms %3i %3i %8.6f \n',h,m,s);
           fprintf(1,'utc %8.6f ',utc );
           [h,m,s] = sec2hms( utc );
           fprintf(1,'hms %3i %3i %8.6f \n',h,m,s);
           fprintf(1,'tai %8.6f',tai );
           [h,m,s] = sec2hms( tai );
           fprintf(1,'hms %3i %3i %8.6f \n',h,m,s);
           fprintf(1,'tt  %8.6f ttt  %16.12f jdtt  %18.11f ',tt,ttt,jdtt+jdttfrac );
           [h,m,s] = sec2hms( tt );
           fprintf(1,'hms %3i %3i %8.6f \n',h,m,s);
           fprintf(1,'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n',tdb,ttdb,jdtdb+jdtdbfrac );

           fprintf(1,'\n start processing, start from the itrf vector: \n');
           fprintf(1,' ecef %14.7f%14.7f%14.7f',recef );
           fprintf(1,' v %14.9f%14.9f%14.9f',vecef );
           fprintf(1,' a %14.9f%14.9f%14.9f\n',aecef );
           
           % -------- ecef2teme    - transform ecef to teme vectors
           [rteme,vteme,ateme] = ecef2teme  ( recef,vecef,aecef,ttt,jdut1+jdut1frac,lod,xp,yp );
           fprintf(1,'ecef-teme %14.7f%14.7f%14.7f',rteme );
           fprintf(1,' v %14.9f%14.9f%14.9f',vteme );
           fprintf(1,' a %14.9f%14.9f%14.9f\n',ateme );

           % -------- teme2ecef    - transform teme to ecef vectors
           [rj2,vj2,aj2] = teme2ecef( rteme,vteme,ateme,ttt,jdut1+jdut1frac,lod,xp,yp );
           fprintf(1,'teme-ecef %14.7f%14.7f%14.7f',rj2 );
           fprintf(1,' v %14.9f%14.9f%14.9f',vj2 );
           fprintf(1,' a %14.9f%14.9f%14.9f\n',aj2 );

clear all;

%    sgp4 example
           ateme = [0.001;0.002;0.003];
           typerun = 'v';
           if (typerun == 'm')
               typeinput = input('input mfe, epoch (YMDHMS), or dayofyr approach, m,e,d: ','s');
           else
                typeinput = 'e';
           end;
           whichconst = 72;
           [tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );

           longstr1 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753';
           longstr2 = '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667     0.00      4320.0        360.00  ';


%           // convert the char string to sgp4 elements
%           // includes initialization of sgp4
           [satrec, startmfe, stopmfe, deltamin] = twoline2rv( whichconst, ...
                        longstr1, longstr2, typerun, typeinput);
           fprintf(1,' %d\n', satrec.satnum);

%           // call the propagator to get the initial state vector value
           [satrec, ro ,vo] = sgp4 (satrec,  0.0);
           fprintf(1, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n',...
                  satrec.t,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));

           tsince = stopmfe;
           tsince = 0.0;
           [satrec, rteme, vteme] = sgp4 (satrec,  tsince);
           fprintf(1, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n',...
                  satrec.t,rteme(1),rteme(2),rteme(3),vteme(1),vteme(2),vteme(3));

           % now convert vectors
           % -------- teme2ecef    - transform teme to ecef vectors
           jdutc = satrec.jdsatepoch;
           jdutcfrac = satrec.jdsatepochf + tsince/1440.0;
           [year,mon,day,hr,min,sec] = invjday(jdutc, jdutcfrac);
           
           % if assume tle epoch is utc and dut1, xp, etc are zero,
%           dut1 = 0.0;
%           dat  = 32;
%           xp   =  0.0;  % "
%           yp   =  0.0;
%           lod  =  0.0;
%           ddpsi = 0.0;  % "
%           ddeps = 0.0;

           % if have values for the date in question, 
           dut1 = 0.2048315;
           dat  = 32;
           xp   =  0.109600;  % "
           yp   =  0.284144;
           lod  =  0.0004116;
           ddpsi = -0.054055;  % "
           ddeps = -0.006183;

           timezone = 0;
           [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
           = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

           [recef,vecef,aecef] = teme2ecef( rteme',vteme',ateme,ttt,jdut1+jdut1frac,lod,xp,yp );
           fprintf(1,'teme-ecef %14.7f%14.7f%14.7f',recef );
           fprintf(1,' v %14.9f%14.9f%14.9f',vecef );
           fprintf(1,' a %14.9f%14.9f%14.9f\n',aecef );

       
           r = recef;
           v= vecef;

           magr = mag(r);
           temp= sqrt( r(1)*r(1) + r(2)*r(2) );
           if ( temp < 0.00000001 )
               lon = atan2( v(2) , v(1) );
             else
               lon = atan2( r(2) , r(1) );
             end
           lat = asin( r(3)/magr );

           alt = magr - radiusearthkm;
           rad = 180.0/pi;
           fprintf(1,' lat %14.9f  lon %14.9f alt %14.9f\n',lat*rad, lon*rad, alt );


