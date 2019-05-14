%     -----------------------------------------------------------------
%
%                              ex3_1415.m
%
%  this file tests the reduction functions.
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
%            30 mar 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************

        % GEO test
%        recef = [24796.91929150; -34115.87092340; 10.22606210];    
%        vecef = [-0.0009791780;	-0.0014765380;	-0.0009287760];
%        aecef = [0.001;0.002;0.003];
        % LEO test
        recef = [-1033.4793830;  7901.2952754;  6380.3565958];    
        vecef = [-3.225636520;  -2.872451450;   5.531924446];
        aecef = [0.001;0.002;0.003];

        conv = pi / (180.0*3600.0);
        
        year=2004;
        mon = 4;
        day = 6;
        hr =  7;
        min= 51;
        sec= 28.386009;

        dut1 = -0.4399619;  % sec
        dat  = 32;         % sec
        xp   = -0.140682 * conv;  % " to rad
        yp   =  0.333309 * conv;
        lod  =  0.0015563;
        ddpsi = -0.052195 * conv;  % " to rad
        ddeps = -0.003875 * conv;
        ddx = -0.000205 * conv;  % " to rad
        ddy = -0.000136 * conv;

        
%         % --------------- 12 convtime test ------------------------
%            recef = [-1060.0995722 ;  -409.3834569  ; 6881.3501455   ];
%            vecef = [  1.74326199 ;  7.38336350 ;  0.70878112 ];
%            aecef = [0.001;0.002;0.003];
% 
%            year=2010;
%            mon = 8;
%            day = 13;
%            hr =  2;
%            min= 52;
%            sec= 39.243;
% %2010 08 09 55417  0.156184  0.457790 -0.0470509  0.0001382 -0.073197 -0.010004 -0.000051 -0.000001  34 
% %2010 08 10 55418  0.158482  0.456716 -0.0473484  0.0004451 -0.072985 -0.009943 -0.000028 -0.000046  34 
% %2010 08 11 55419  0.161130  0.455342 -0.0479526  0.0006765 -0.073110 -0.009650  0.000001  0.000007  34 
% 
%            dut1 = -0.0473484;
%            dat  = 34;
%            xp   = 0.158482 * conv;  % rad
%            yp   =  0.456716 * conv;
%            lod  =  0.0004451;
%            ddpsi = -0.072985 * conv;  % rad
%            ddeps = -0.009943 * conv;
        
        
        timezone=0;
        order = 106;
        eqeterms = 2; % use the extra eqeq terms in j2000
        opt = 'a'; % specify the iau00 approach

        fprintf(1, 'test program for reduction functions \n\n');

        fprintf(1, 'input data \n\n');
        fprintf(1, ' year %5i ', year);
        fprintf(1, ' mon %4i ', mon);
        fprintf(1, ' day %3i ', day);
        fprintf(1, ' %3i:%2i:%8.6f\n ', hr, min, sec );
        fprintf(1, ' dut1 %8.6f s', dut1);
        fprintf(1, ' dat %3i s', dat);
        fprintf(1, ' xp %8.6f "', xp / conv);
        fprintf(1, ' yp %8.6f "', yp / conv);
        fprintf(1, ' lod %8.6f s\n', lod);
        fprintf(1, ' ddpsi %8.6f " ddeps  %8.6f\n', ddpsi/conv, ddeps/conv);
        fprintf(1, ' ddx   %8.6f " ddy    %8.6f\n', ddx/conv, ddy/conv);
        fprintf(1, ' order %3i  eqeterms %31  opt %3s \n', order, eqeterms, opt );
        fprintf(1, 'units are km and km/s and km/s2\n' );

        timezone=0;

        % -------- convtime    - convert time from utc to all the others
        fprintf(1, 'convtime results\n');
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] ...
                    = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        fprintf(1, 'ut1 %8.6f tut1 %16.12f jdut1 %18.11f ', ut1, tut1, jdut1+jdut1frac );
        [h, m, s] = sec2hms( ut1 );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'utc %8.6f ', utc );
        [h, m, s] = sec2hms( utc );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tai %8.6f', tai );
        [h, m, s] = sec2hms( tai );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tt  %8.6f ttt  %16.12f jdtt  %18.11f ', tt, ttt, jdtt+jdttfrac );
        [h, m, s] = sec2hms( tt );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n', tdb, ttdb, jdtdb+jdtdbfrac );

        % -------- precess     - transformation matrix for precession
        fprintf(1, 'precession matrix \n');
        [prec, psia, wa, ea, xa] = precess ( ttt, '80' );
        fprintf(1, '%16.11f %16.11f %16.11f\n', prec );

        % -------- nutation    - transformation matrix for nutation
        fprintf(1, 'nutation matrix \n');
        [deltapsi00, trueeps00, meaneps00, omega00, nut00] = nutation(ttt, 0, 0);
        fprintf(1, '%16.11f %16.11f %16.11f\n', nut00 );
        [deltapsi, trueeps, meaneps, omega, nut] = nutation(ttt, ddpsi, ddeps);
        fprintf(1, '%16.11f %16.11f %16.11f\n', nut );

        % -------- sidereal    - transformation matrix for sidereal time
        fprintf(1, 'sidereal time matrix \n');
        [st00, stdot00] = sidereal(jdut1+jdut1frac, deltapsi, meaneps, omega, lod, 0);
        fprintf(1, '%16.11f %16.11f %16.11f\n', st00 );
        [st, stdot] = sidereal(jdut1+jdut1frac, deltapsi, meaneps, omega, lod, 2);
        fprintf(1, '%16.11f %16.11f %16.11f\n', st );
% note you could also have two calculations with the wo corr nutation values



        % -------- polarm      - transformation matrix for polar motion
        fprintf(1, 'polar motion matrix \n');
        [pm] = polarm(xp, yp, ttt, '80');
        fprintf(1, '%16.11f %16.11f %16.11f\n', pm );

        %-------------------
        % sample runs for book
        [recigg, vecigg, aecig] = iau00f2iS ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f', recigg );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecigg );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecig );

        [recig, vecig, aecig] = ecef2eci(recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps);
        fprintf(1, 'GCRF 2 w corr IAU-76/FK5   %14.7f %14.7f %14.7f', recig );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecig );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecig );
        %-------------------
        
  
        
        fprintf(1, '\n\n ============== convert various coordinate systems from ecef =================== \n');
        fprintf(1, 'ITRF          IAU-76/FK5   %14.7f %14.7f %14.7f', recef );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecef );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecef );

        % -------- pef transformations
        [rpef, vpef, apef] = ecef2pef  ( recef, vecef, aecef, xp, yp, ttt );
        fprintf(1, 'PEF           IAU-76/FK5   %14.7f %14.7f %14.7f', rpef );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vpef );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', apef );
        
        % -------- tod transformations
        [rtod, vtod, atod] = ecef2tod  ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps );
        fprintf(1, 'TOD 2 w corr  IAU-76/FK5   %14.7f %14.7f %14.7f', rtod );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vtod );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', atod );
        [rtod20, vtod20, atod20] = ecef2tod  ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 2, 0, 0 );
        fprintf(1, 'TOD 2 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f', rtod20 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vtod20 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', atod20 );
        [rtod00, vtod00, atod00] = ecef2tod  ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 0, 0, 0 );
        fprintf(1, 'TOD 0 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f', rtod00 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vtod00 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', atod00 );
        
        % -------- mod transformations
        [rmod, vmod, amod] = ecef2mod  ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps );
        fprintf(1, 'MOD 2 w corr  IAU-76/FK5   %14.7f %14.7f %14.7f', rmod );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vmod );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', amod );
        [rmod20, vmod20, amod20] = ecef2mod  ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 2, 0, 0 );
        fprintf(1, 'MOD 2 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f', rmod20 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vmod20 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', amod20 );
        [rmod00, vmod00, amod00] = ecef2mod  ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 0, 0, 0 );
        fprintf(1, 'MOD 0 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f', rmod00 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vmod00 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', amod00 );

        % -------- eci transformations
        [recig, vecig, aecig] = ecef2eci(recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps);
        fprintf(1, 'GCRF 2 w corr IAU-76/FK5   %14.7f %14.7f %14.7f', recig );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecig );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecig );
        % now do check back to itrf to make sure you get the same vector
        [recef2w, vecef2w, aecef2w] = eci2ecef(recig, vecig, aecig, ttt, jdut1, lod, xp, yp, 2, ddpsi, ddeps);
        fprintf(1, 'eci-ecef 2, w corr         %14.7f %14.7f %14.7f', recef2w );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecef2w );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecef2w );

        [reci20, veci20, aeci20] = ecef2eci(recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 2, 0, 0);
        fprintf(1, 'GCRF 2wo corr IAU-76/FK5   %14.7f %14.7f %14.7f', reci20 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', veci20 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aeci20 );
        [recef, vecef, aecef] = eci2ecef(reci20, veci20, aeci20, ttt, jdut1+jdut1frac, lod, xp, yp, 2, 0, 0);
        fprintf(1, 'ecef 2wo corr              %14.7f %14.7f %14.7f', recef );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecef );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecef );
        
        [reci00, veci00, aeci00] = ecef2eci(recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 0, 0, 0);
        fprintf(1, 'GCRF 0wo corr IAU-76/FK5   %14.7f %14.7f %14.7f', reci00 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', veci00 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aeci00 );
        [recef, vecef, aecef] = eci2ecef(reci00, veci00, aeci00, ttt, jdut1+jdut1frac, lod, xp, yp, 0, 0, 0);
        fprintf(1, 'ecef 0 wo corr             %14.7f %14.7f %14.7f', recef );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecef );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecef );

        
        
        % -------- ecef2teme    - transform ecef to teme vectors
        [rteme, vteme, ateme] = ecef2teme  ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp );
        fprintf(1, 'TEME                       %14.7f %14.7f %14.7f', rteme );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vteme );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', ateme );
        [receft, veceft, aeceft] = teme2ecef  ( rteme, vteme, ateme, ttt, jdut1+jdut1frac, lod, xp, yp );
        fprintf(1, 'teme - ecef                %14.7f %14.7f %14.7f', receft );
        fprintf(1, ' v %14.9f %14.9f %14.9f', veceft );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aeceft );


        tm = fk4;
        r1950 = tm' * reci20;
        v1950 = tm' * veci20;
        fprintf(1, 'eci - FK4                  %14.7f %14.7f %14.7f', r1950 );
        fprintf(1, ' v %14.9f %14.9f %14.9f \n', v1950 );

        reci20i = tm * r1950;
        veci20i = tm * v1950;
        fprintf(1, 'FK4 - eci                  %14.7f %14.7f %14.7f', reci20i );
        fprintf(1, ' v %14.9f %14.9f %14.9f \n', veci20i );
        
        % -------- pef2eci     - transform pef to eci vectors
        [reci, veci, aeci] = pef2eci(rpef, vpef, apef, ttt, jdut1+jdut1frac, lod, 2, 0, 0);
        fprintf(1, 'pef-eci 20                 %14.7f %14.7f %14.7f', reci );
        fprintf(1, ' v %14.9f %14.9f %14.9f', veci );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aeci );

        % -------- tod2eci     - transform tod to eci vectors
        [reci0, veci0, aeci0] = tod2eci(rtod00, vtod00, atod00, ttt, 0, 0);
        fprintf(1, 'tod-eci wo corr            %14.7f %14.7f %14.7f', reci0 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', veci0 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aeci0 );

        % -------- mod2eci     - transform mod to eci vectors
        [recig, vecig, aecig] = mod2eci(rmod, vmod, amod, ttt);
        fprintf(1, 'mod-eci                    %14.7f %14.7f %14.7f', recig );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecig );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecig );

        
        [rtod0, vtod0, atod0] = eci2tod(reci, veci, aeci, ttt, 0, 0);
        fprintf(1, 'eci-TOD wo corr            %14.7f %14.7f %14.7f', rtod0 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vtod0 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', atod0 );

        [rpef20, vpef20, apef20] = eci2pef(reci20, veci20, aeci20, ttt, jdut1, lod, 2, 0, 0);
        fprintf(1, 'eci-PEF 20,                %14.7f %14.7f %14.7f', rpef20 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vpef20 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', apef20 );

        [rpef00, vpef00, apef00] = eci2pef(reci00, veci00, aeci00, ttt, jdut1+jdut1frac, lod, 0, 0, 0);
        fprintf(1, 'eci-PEF 00                 %14.7f %14.7f %14.7f', rpef00 );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vpef00 );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', apef00 );

        [rtodwc, vtodwc, atodwc] = eci2tod(recig, vecig, aecig, ttt, ddpsi, ddeps);
        fprintf(1, 'eci-TOD wc                 %14.7f %14.7f %14.7f', rtodwc );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vtodwc );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', atodwc );

        [rpef, vpef, apef] = eci2pef(recig, vecig, aecig, ttt, jdut1+jdut1frac, lod, 2, ddpsi, ddeps);
        fprintf(1, 'eci-pef 2 wc               %14.7f %14.7f %14.7f', rpef );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vpef );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', apef );



%        % -------- eci2teme    - transform eci to teme vectors
%        [rteme, vteme, ateme] = eci2teme(reci, veci, aeci, ttt, 4, 0, opt);
%        fprintf(1, 'order   4  terms   0  opt %3s \n', opt );
%        fprintf(1, 'eci-teme %14.7f %14.7f %14.7f', rteme );
%        fprintf(1, ' v %14.9f %14.9f %14.9f', vteme );
%        fprintf(1, ' a %14.9f %14.9f %14.9f\n', ateme );



        % --------------------------- now do iau2000
        % -------- eci2ecef    - transform iau2000 eci to ecef vectors
%    recef=[ -1033.4793830  7901.2952758  6380.3565953];
%    vecef=[    -3.225636520  -2.872451450   5.531924446];
        fprintf(1, '\n\n\n\n now do iau2000 processing., start from the itrf vector: \n');
        fprintf(1, 'eci-ecef                   %14.7f %14.7f %14.7f', recef );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecef );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecef );

        % -------- ecef2eci    - transform iau2000 ecef to eci vectors
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'a', 0, 0 );
        fprintf(1, 'GCRF          IAU-2006     %14.7f %14.7f %14.7f', recigg );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecigg );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecig );

        [recefg, vecefg, aecefg] = iau00i2f  ( recigg, vecigg, aecig, ttt, jdut1+jdut1frac, lod, xp, yp, 'a', 0, 0 );
        fprintf(1, 'itrf                       %14.7f %14.7f %14.7f', recefg );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecefg );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecefg );

        % -------- ecef2eci    - transform iau2000 ecef to eci vectors
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'b', 0, 0 );
        fprintf(1, 'GCRF          IAU-2000 B   %14.7f %14.7f %14.7f', recigg );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecigg );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecig );

        [recefg, vecefg, aecefg] = iau00i2f ( recigg, vecigg, aecig, ttt, jdut1+jdut1frac, lod, xp, yp, 'b', 0, 0  );
        fprintf(1, 'itrf                       %14.7f %14.7f %14.7f', recefg );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecefg );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecefg );

        % -------- ecef2eci    - transform iau2000 ecef to eci vectors
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f', recigg );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecigg );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecig );

        [recefg, vecefg, aecefg] = iau00i2f  ( recigg, vecigg, aeci, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
        fprintf(1, 'itrf                       %14.7f %14.7f %14.7f', recefg );
        fprintf(1, ' v %14.9f %14.9f %14.9f', vecefg );
        fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecefg );


        fprintf(1, 'differences \n\n');
        fprintf(1, ' gcrf-fk5  %14.9fm \n', mag(recig-reci)*1000 );
        fprintf(1, ' fk5-mod   %14.9fm \n', mag(reci-rmod)*1000 );
        fprintf(1, ' mod-tod   %14.9fm \n', mag(rmod-rtod)*1000 );
        fprintf(1, ' tod-pef   %14.9fm \n', mag(rtod-rpef)*1000 );
        fprintf(1, ' pef-ecef  %14.9fm \n', mag(rpef-recef)*1000 );
        fprintf(1, ' ggcrf-fk5 %14.9fm \n', mag(recigg-reci)*1000 );

        
        
        % do sensitivity tests        
        % GEO tests -----------------------------------------------------------
        recef = [24796.91929150; -34115.87092340; 10.22606210];    
        vecef = [-0.0009791780;	-0.0014765380;	-0.0009287760];
        aecef = [0.001;0.002;0.003];
        ddxerror = 0.0001;

        % find answer one time
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigans, vecigans, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
        fprintf(1, 'GCRF ans      IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigans, vecigans );
        
        % ---- DAT err of 1 sec
        dat = 31;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c1 dat 1 s GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        dat = 32;

        % ---- TT off by a min
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        dat = dat - 60.0;
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c2 tt 1 m GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        dat = dat + 60.0;
       
        % ---- dpsi corretion off by 0.001 arcsec      
        ddx = ddx - ddxerror * conv;  % " to rad
        ddy = ddy - ddxerror * conv;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c6 ddxy GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        ddx = ddx + ddxerror * conv;  % " to rad
        ddy = ddy + ddxerror * conv;

        % ---- DUT1 off by 0.01 sec
        dut1 = dut1-0.01;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c3 dut1 0.01s GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        dut1 = dut1+0.01;

        % ---- DUT1 off by 0.25 sec
        dut1 = dut1-0.25;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c4 dut1 0.25s GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        dut1 = dut1+0.25;

        % ---- UTC error by 1 sec
        sec = sec - 1;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c5 utc 1s GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        sec = sec + 1;
        fprintf(1, 'ut1 %8.6f tut1 %16.12f jdut1 %18.11f ', ut1, tut1, jdut1+jdut1frac );
        [h, m, s] = sec2hms( ut1 );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'utc %8.6f ', utc );
        [h, m, s] = sec2hms( utc );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tai %8.6f', tai );
        [h, m, s] = sec2hms( tai );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tt  %8.6f ttt  %16.12f jdtt  %18.11f ', tt, ttt, jdtt+jdttfrac );
        [h, m, s] = sec2hms( tt );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n', tdb, ttdb, jdtdb+jdtdbfrac );

        
        % LEO tests -----------------------------------------------------------
        recef = [-1033.4793830;  7901.2952754;  6380.3565958;];    
        vecef = [-3.225636520;  -2.872451450;   5.531924446;];
        aecef = [0.001;0.002;0.003];
        
%           recef = [3961.00354581;  6010.75117182; 4619.30093637 ];
%           vecef = [-5.315109071309; 3.96813865785; 1.752758565833 ];

        % find answer one time
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigans, vecigans, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
        fprintf(1, '\n\nGCRF ans      IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigans, vecigans );
        
        % ---- DAT err of 1 sec
        dat = 31;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c1 dat 1 s GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        dat = 32;

        % ---- TT off by a min
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        dat = dat - 60.0;
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c2 tt 1 m GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        dat = dat + 60.0;
       
        % ---- dpsi corretion off by 0.001 arcsec      
        ddx = ddx - ddxerror * conv;  % " to rad
        ddy = ddy - ddxerror * conv;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c6 ddxy GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        ddx = ddx + ddxerror * conv;  % " to rad
        ddy = ddy + ddxerror * conv;

        % ---- DUT1 off by 0.01 sec
        dut1 = dut1-0.01;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c3 dut1 0.01s GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        dut1 = dut1+0.01;

        % ---- DUT1 off by 0.25 sec
        dut1 = dut1-0.25;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c4 dut1 0.25s GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        dut1 = dut1+0.25;

        % ---- UTC error by 1 sec
        sec = sec - 1;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [recigg, vecigg, aecig] = iau00f2i ( recef, vecef, aecef, ttt, jdut1, lod, xp, yp, 'c', ddx, ddy );
%        fprintf(1, 'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f  v %14.9f %14.9f %14.9f\n', recigg, vecigg );
        fprintf(1, ' c5 utc 1s GCRF %14.9fm \n', mag(recigans-recigg)*1000 );
        sec = sec + 1;
        fprintf(1, 'ut1 %8.6f tut1 %16.12f jdut1 %18.11f ', ut1, tut1, jdut1+jdut1frac );
        [h, m, s] = sec2hms( ut1 );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'utc %8.6f ', utc );
        [h, m, s] = sec2hms( utc );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tai %8.6f', tai );
        [h, m, s] = sec2hms( tai );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tt  %8.6f ttt  %16.12f jdtt  %18.11f ', tt, ttt, jdtt+jdttfrac );
        [h, m, s] = sec2hms( tt );
        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
        fprintf(1, 'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n', tdb, ttdb, jdtdb+jdtdbfrac );
        
        
        
