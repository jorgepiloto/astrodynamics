        %
        %  exsgp4_teme
        %
        %  tests the teme conversions including of date and of epoch
        %

        conv = pi / (180*3600);

        % ---------- test the teme to pef conversion ------------
        recef = [-1033.4793830;  7901.2952754;  6380.3565958];    
        vecef = [-3.225636520;  -2.872451450;   5.531924446];
        aecef = [1; 1; 1];
       
        year = 2004;
        mon  =   4;
        day  =   6;
        hr   =   7;
        min  =  51;
        sec  =  28.386009;
        dut1 = -0.4399619;
        dat  = 32;
        xp   = -0.140682 * conv;  % " to rad
        yp   =  0.333309 * conv;
        lod  =  0.0015563;  % sec
        ddpsi = -0.052195 * conv;  % " to rad
        ddeps = -0.003875 * conv;
        ddx = -0.000205 * conv;  % " to rad
        ddy = -0.000136 * conv;
        timezone=0;
        order = 106;
        eqeterms = 2; % use the kinematic eqe terms after 1997
        
        fprintf(1,'test program for reduction functions \n\n');
        fprintf(1,'input data \n\n');
        fprintf(1,'year %5i ',year);
        fprintf(1,' mon %4i ',mon);
        fprintf(1,' day %3i ',day);
        fprintf(1,'%3i:%2i:%8.6f\n ',hr,min,sec );
        fprintf(1,'dut1 %8.6f s',dut1);
        fprintf(1,' dat %3i s',dat);
        fprintf(1,' xp %8.6f "',xp/conv);
        fprintf(1,' yp %8.6f "',yp/conv);
        fprintf(1,' lod %8.6f s\n',lod);
        fprintf(1,' ddpsi %8.6f " ddeps  %8.6f\n',ddpsi/conv,ddeps/conv);
        fprintf(1,' ddx %8.6f " ddy  %8.6f\n',ddx/conv,ddy/conv);
        fprintf(1,'order %3i  eqeterms %31  \n',order,eqeterms );
        fprintf(1,'units are km and km/s \n' );

        timezone=0;

        % -------- convtime    - convert time from utc to all the others
        fprintf(1,'convtime results\n');
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
              = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        fprintf(1,'ut1 %8.6f tut1 %16.12f jdut1 %18.11f \n',ut1,tut1,jdut1 );
        
        % ---------------------- teme transformations ---------------------
        % -------- ecef2teme    - transform ecef to teme vectors
        [rteme, vteme, ateme] = ecef2teme(recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, eqeterms);
        fprintf(1,'\n\n start from ecef \n');
        fprintf(1,'ecef-teme\n rteme %14.7f %14.7f %14.7f',rteme );
        fprintf(1,' vteme %14.9f %14.9f %14.9f\n',vteme );

        % -------- teme2ecef    - transform teme to ecef vectors
        [recef1, vecef1, aecef1] = teme2ecef(rteme, vteme, ateme, ttt, jdut1+jdut1frac, lod, xp, yp, eqeterms);
        fprintf(1,'teme-ecef\n recef %14.7f %14.7f %14.7f',recef1 );
        fprintf(1,' vecef %14.9f %14.9f %14.9f\n',vecef1 );
        dr = 1000*(recef - recef1);  % in m
        fprintf(1,'diff in teme %14.7f %14.7f %14.7f %14.7f \n',dr(1),dr(2),dr(3),mag(dr) );

        % ---------------------- eci transformations ----------------------
        % -------- teme2eci    - transform teme to eci vectors
        fprintf(1,'\n\n start from teme \n');
        [reci, veci, aeci] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);
        fprintf(1,'teme-eci\n reci %14.7f %14.7f %14.7f',reci );
        fprintf(1,' veci %14.9f %14.9f %14.9f\n',veci );

        % -------- eci2teme    - transform eci to teme vectors
        [rteme1, vteme1, ateme1] = eci2teme(reci, veci, aeci, ttt, ddpsi, ddeps);
        fprintf(1,'ecef-teme\n rteme %14.7f %14.7f %14.7f',rteme1 );
        fprintf(1,' vteme %14.9f %14.9f %14.9f\n',vteme1 );
        dr = 1000*(rteme - rteme1);  % in m
        fprintf(1,'diff in teme %14.7f %14.7f %14.7f %14.7f \n',dr(1),dr(2),dr(3),mag(dr) );
        
        % check standard ecef to eci transformation
        [recig, vecig, aecig] = ecef2eci(recef, vecef, aecef, ttt, jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps);
        fprintf(1, 'GCRF 2 w corr IAU-76/FK5   %14.7f %14.7f %14.7f', recig );
        fprintf(1, ' v %14.9f %14.9f %14.9f\n', vecig );
        dr = 1000*(reci - recig);  % in m
        fprintf(1,'diff in eci %14.7f %14.7f %14.7f %14.7f \n',dr(1),dr(2),dr(3),mag(dr) );

        pause;      

        % -------------------- teme of epoch/date --------------------
        fprintf(1,'\n\n =================== now do teme of date and of epoch example  ================= \n');        
        %       //typerun = 'c' compare 1 year of full satcat data
        %       //typerun = 'v' verification run, requires modified elm file with
        %       //              start stop and delta times
        rad = 180.0/pi;
        opsmode= 'a';  % afspc approach
        whichconst = 72;
        eqeterms = 2;

        longstr1 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753';
        longstr2 = '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667     0.0      4320.0       720.00000';

        %       // convert the char string to sgp4 elements
        %       // includes initialization of sgp4
        typerun = 'v';
        [startmfe, stopmfe, deltamin, satrec] = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);
        fprintf(1,'%11.7f %11.7f %11.7f \n',startmfe, stopmfe, deltamin );           
        fprintf(1,' %d\n', satrec.satnum);

        %      // call the propagator to get the initial state vector value
        [satrec,ro,vo] = sgp4 (satrec,  0.0);

        fprintf(1, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n',...
                satrec.t,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));

        tsince = startmfe;

        %  // check so the first value isn't written twice
        if ( abs(tsince) > 1.0e-8 )
              tsince = tsince - deltamin;
        end
        % eop data
        % date      MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
        %2000 06 27 51722  0.109021  0.285481  0.2052872  0.0004465 -0.053678 -0.006320 -0.000025  0.000028  32
        [year,mon,day,hr,min,sec] = invjday (satrec.jdsatepoch, satrec.jdsatepochf);
       
        dut1 = 0.2052872;
        dat  = 32;
        xp   = 0.109021 * conv;  % " to rad
        yp   = 0.285481 * conv;
        lod  =  0.0004465;  % sec
        ddpsi = -0.053678 * conv;  % " to rad
        ddeps = -0.006320 * conv;
        ddx = -0.000025 * conv;  % " to rad
        ddy = 0.000028 * conv;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
              = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

        fprintf(1,'year %5i ',year);
        fprintf(1,' mon %4i ',mon);
        fprintf(1,' day %3i ',day);
        fprintf(1,'%3i:%2i:%8.6f\n',hr,min,sec );
          
         % for teme of epoch, find the transformation matrix at t0
        [prec,psia,wa,ea,xa] = precess ( ttt, '80' );
        [deltapsie, trueepse, meanepse, omegae, nute] = nutation(ttt, ddpsi, ddeps);  % , ddpsi, ddeps);

        eqeg = deltapsie * cos(meanepse);
        eqeg = rem (eqeg, 2.0*pi);

        eqe(1,1) =  cos(eqeg);
        eqe(1,2) =  sin(eqeg);
        eqe(1,3) =  0.0;
        eqe(2,1) = -sin(eqeg);
        eqe(2,2) =  cos(eqeg);
        eqe(2,3) =  0.0;
        eqe(3,1) =  0.0;
        eqe(3,2) =  0.0;
        eqe(3,3) =  1.0;

        tme = prec * nute * eqe';
        
        
        % // loop to perform the propagation
        % while ((tsince < stopmfe) && (satrec.error == 0))
        while ((tsince < stopmfe))
            tsince = tsince + deltamin;

            if(tsince > stopmfe)
                tsince = stopmfe;
            end

            [satrec,ro,vo] = sgp4 (satrec,  tsince);

            if (satrec.error == 0)
                if ((typerun ~= 'e') && (typerun ~= 'd'))
                   jd = satrec.jdsatepoch;
                   jdf = satrec.jdsatepochf + tsince/1440.0;
                   [year,mon,day,hr,min,sec] = invjday ( jd, jdf );

                   fprintf(1, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f  %5i%3i%3i %2i:%2i:%9.6f\n',...
                           satrec.t,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3), year,mon,day,hr,min,sec);
                 else
                   fprintf(1, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f',...
                           tsince,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));
                   [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo);
                   fprintf(1, ' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f \n',...
                              a, ecc, incl*rad, omega*rad, argp*rad, nu*rad, m*rad);
                 end
            end %// if satrec.error == 0

        end %// while propagating the orbit

        rteme = [ro(1); ro(2); ro(3)];       
        vteme = [vo(1); vo(2); vo(3)];       
        fprintf(1,' rteme %14.7f %14.7f %14.7f',rteme );
        fprintf(1,' vteme %14.9f %14.9f %14.9f\n',vteme );

        % for teme of date, find the transformation matrix at t=t0 + 3days
        %2000 06 30 51725  0.110329  0.281805  0.2042651  0.0001678 -0.054522 -0.006209  0.000002  0.000016  32
        dut1 = 0.2042651;
        dat  = 32;
        xp   = 0.110329 * conv;  % " to rad
        yp   = 0.281805 * conv;
        lod  =  0.0001678;  % sec
        ddpsi = -0.054522 * conv;  % " to rad
        ddeps = -0.006209 * conv;
        ddx = 0.000002 * conv;  % " to rad
        ddy = 0.000016 * conv;
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
              = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        
        % ------ teme of date
        [prec,psia,wa,ea,xa] = precess ( ttt, '80' );
        [deltapsi, trueeps, meaneps, omega, nut] = nutation(ttt, ddpsi, ddeps);  % , ddpsi, ddeps);

        eqeg = deltapsi * cos(meaneps);
        eqeg = rem (eqeg, 2.0*pi);

        eqe(1,1) =  cos(eqeg);
        eqe(1,2) =  sin(eqeg);
        eqe(1,3) =  0.0;
        eqe(2,1) = -sin(eqeg);
        eqe(2,2) =  cos(eqeg);
        eqe(2,3) =  0.0;
        eqe(3,1) =  0.0;
        eqe(3,2) =  0.0;
        eqe(3,3) =  1.0;

        tm = prec * nut * eqe';
        fprintf(1,'of date transformation teme - eci   %11.7f  %11.7f  %11.7f \n', tm);
        % -------- teme2ecef    - transform teme to ecef vectors
        [recefd, vecefd, aecefd] = teme2ecef(rteme, vteme, ateme, ttt, jdut1+jdut1frac, lod, xp, yp, eqeterms);
        fprintf(1,'teme-ecef of date\n recefd %14.7f %14.7f %14.7f',recefd );
        fprintf(1,' vecefd %14.9f %14.9f %14.9f\n',vecefd );
          
        % -------- teme2eci    - transform teme to eci vectors
        fprintf(1,'\n\n start from teme \n');
        [recid, vecid, aecid] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);
        fprintf(1,'teme-eci of date\n recid %14.7f %14.7f %14.7f',recid );
        fprintf(1,' vecid %14.9f %14.9f %14.9f\n',vecid );

        % convert ecef to eci standard
        % -------- ecef2eci    - transform ecef to eci vectors
        [reci, veci, aeci] = ecef2eci(recefd, vecefd, aecefd, ttt, jdut1+jdut1frac, lod, xp, yp, 2,  ddpsi, ddeps);
        fprintf(1,'ecef-eci fk5/iau76\n reci %14.7f%14.7f%14.7f',reci );
        fprintf(1,' veci %14.9f%14.9f%14.9f\n',veci );


        % ------ teme of epoch
        % -------- teme2eci    - transform teme to eci vectors using epoch tm matrix
        fprintf(1,'of epoch transformation teme - eci   %11.7f  %11.7f  %11.7f \n', tme);

        recie = tme * rteme;
        vecie = tme * vteme;
        aecie = tme * ateme;

        fprintf(1,'teme-eci of epoch\n recie %14.7f %14.7f %14.7f',recie );
        fprintf(1,' vecie %14.9f %14.9f %14.9f\n',vecie );
        
        % ---- now convert back to ecef using std techniques
        [recefe,vecefe,aecefe] = eci2ecef( recie, vecie, aecie, ttt, jdut1+jdut1frac, lod, xp, yp, eqeterms, ddpsi, ddeps );
                
        fprintf(1,'teme-ecef of epoch\n recefe %14.7f %14.7f %14.7f',recefe );
        fprintf(1,' vecefe %14.9f %14.9f %14.9f\n',vecefe );
       
        dr = recid - recie;
        dv = vecid - vecie;
        fprintf(1,'eci diffs\n dr (m) %14.9f %14.9f %14.9f',dr*1000.0 );
        fprintf(1,' dv %14.9f %14.9f %14.9f\n',dv*1000.0 );
        fprintf(1,'%11.7f \n',mag(dr));
     
        dr = recefd - recefe;
        dv = vecefd - vecefe;
        fprintf(1,'ecef diffs\n dr (m) %14.9f %14.9f %14.9f',dr*1000.0 );
        fprintf(1,' dv %14.9f %14.9f %14.9f\n',dv*1000.0 );
        fprintf(1,'%11.7f \n',mag(dr));
        
        
        
        
        
  