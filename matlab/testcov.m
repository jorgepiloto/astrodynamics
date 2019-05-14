%
%     -----------------------------------------------------------------
%
%                               testcov.m
%
%  this file tests the accuracy of the covariance functions.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                2013
%                          by david vallado
%
%     (h)               email davallado@gmail.com
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            22 jun 15  david vallado
%                         many fixes for new paper
%  changes :
%             9 aug 03  david vallado
%                         fix units on n in equinoc?
%            25 jul 03  david vallado
%                         fixes for tm
%            20 jul 03  david vallado
%                         fixes to add extras for paper
%            25 jun 03  david vallado
%                         update for alternate cases
%            26 may 03  david vallado
%                         fix units and values
%            28 oct 02  david vallado
%                         fix covariance transformations
%             5 sep 02  david vallado
%                         fix cov trace input, misc
%            26 aug 02  david vallado
%                         work on partials for covariance and test cases
%            17 aug 02  david vallado
%                         work on partials for covariance
%            30 jun 02  david vallado
%                         breakout reduction, time, 2body
%            24 may 02  david vallado
%                         original baseline
%
%     *****************************************************************

        constastro;

        small = 1.0e-18;
        rad   = 180.0/pi;
        rad2  = rad*rad;

        anom = 'true';
        anom = 'mean';
        testnum = 3;

        ddpsi = 0.0;
        ddeps = 0.0;

        % --------------- 0. AF test
        if testnum==0
           reci = [3961.74426025;  6010.21561092; 4619.36257583 ];
           veci = [-5.314643385545; 3.964357584990; 1.752939152769 ];
           aeci = [0.001;0.002;0.003];

           year = 2000;
           mon  =   6;
           day  =  28;
           hr   =  15;
           min  =   8;
           sec  =  51.655;

           dut1 =  0.16236;
           dat  = 21;
           xp   =  0.0987;
           yp   =  0.2860;
           lod  =  0.0;
           timezone = 0;
           order = 4;
           terms = 2;
        end;

        % ------- navy test
        if testnum == 1
           reci = [-605.79221660; -5870.22951108; 3493.05319896;];
           veci = [-1.56825429; -3.70234891; -6.47948395; ];
           aeci = [0.001;0.002;0.003];
           year = 2000;
           mon  =  12;
           day  =  15;
           hr   =  16;
           min  =  58;
           sec  =  50.208;
           dut1 =  0.10597;
           dat  = 32;
           xp   =  0.0;
           yp   =  0.0;
           lod  =  0.0;
           terms = 2;
           timezone= 0;
           order = 106;
           cartcov = [ ...
                       8.04204e-13,-7.41923e-13, 2.02623e-12,-3.78793e-16, 1.77302e-13,-2.48352e-13; ...
                      -7.41923e-13, 2.19044e-12, 3.82034e-12,-1.19010e-15, 2.83844e-13,-3.67925e-13; ...
                       2.02623e-12, 3.82034e-12, 1.97570e-11,-9.35438e-15,-1.47386e-12, 3.01542e-12; ...
                      -3.78793e-16,-1.19010e-15,-9.35438e-15, 1.56236e-17, 8.86093e-17, 4.56974e-16; ...
                       1.77302e-13, 2.83844e-13,-1.47386e-12, 8.86093e-17, 9.67797e-13, 7.23072e-13; ...
                      -2.48352e-13,-3.67925e-13, 3.01542e-12, 4.56974e-16, 7.23072e-13, 1.84123e-12];
                    %  3.06012e-12  1.16922e-11  6.72154e-11  1.28647e-13  4.11455e-13 -3.89727e-12  2.20508e-05
        end

        % ------- accuracy sensitivity test
        if testnum == 2
           reci = [-605.79221660; -5870.22951108; 3493.05319896;];
           veci = [-1.56825429; -3.70234891; -6.47948395; ];
           aeci = [0.001;0.002;0.003];

           year = 2000;
           mon  =  12;
           day  =  15;
           hr   =  16;
           min  =  58;
           sec  =  50.208;
           dut1 =  0.10597;
           dat  = 32;
           xp   =  0.0;
           yp   =  0.0;
           lod  =  0.0;
           timezone= 0;
           terms = 2;
           order = 106;
           covopntw = [ ...
                       1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
                       0.0, 10.0, 0.0, 0.0, 0.0, 0.0; ...
                       0.0, 0.0, 1.0, 0.0, 0.0, 0.0; ...
                       0.0, 0.0, 0.0, 1.0e-6, 0.0, 0.0; ...
                       0.0, 0.0, 0.0, 0.0, 1.0e-4, 0.0; ...
                       0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-6];
           [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
           = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

           [cartstate,classstate,flstate,eqstate] = setcov(reci,veci, ...
            year,mon,day,hr,min,sec,dut1,dat,ttt,jdut1,lod,xp,yp,terms,'y',anom,ddpsi,ddeps);

           printcov( covntw,'ct','m',anom );

           [cartcov,tm] = covo22ct( covopntw,cartstate );
           printcov( cartcov,'ct','m',anom );
pause;
        end

        % ------- trace test
        if testnum == 3
           reci = [-605.79221660; -5870.22951108; 3493.05319896;];
           veci = [-1.56825429; -3.70234891; -6.47948395; ];
           aeci = [0.001;0.002;0.003];
         
           year = 2000;
           mon  =  12;
           day  =  15;
           hr   =  16;
           min  =  58;
           sec  =  50.208;
           dut1 =  0.10597;
           dat  = 32;
           xp   =  0.0;
           yp   =  0.0;
           lod  =  0.0;
           timezone= 0;
           terms = 2;
           order = 106;
           cartcov = [ ...
                       0.81,   1.0e-2, 1.0e-2, 1.0e-3, 1.0e-3, 1.0e-3; ...
                       1.0e-2, 0.81,   1.0e-2, 1.0e-3, 1.0e-3, 1.0e-3; ...
                       1.0e-2, 1.0e-2, 0.81,   1.0e-3, 1.0e-3, 1.0e-3; ...
                       1.0e-4, 1.0e-4, 1.0e-4, 0.000001, 1.0e-6, 1.0e-6; ...
                       1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6, 0.000001, 1.0e-6; ...
                       1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6, 1.0e-6, 0.000001];
           cartcov = [ ...
                       1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
                       0.0, 1.0, 0.0, 0.0, 0.0, 0.0; ...
                       0.0, 0.0, 1.0, 0.0, 0.0, 0.0; ...
                       0.0, 0.0, 0.0, 0.000001, 0.0, 0.0; ...
                       0.0, 0.0, 0.0, 0.0, 0.000001, 0.0; ...
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.000001];

           classcovtrace = [ ...
             0.6489403e+01,   0.2979136e-06,  0.2117582e-21,  0.5293956e-21,  0.3465237e-01,  0.9048777e-02; ...
             0.2979136e-06,   0.4353986e-13,  0.6708473e-26,  0.1091152e-25,  0.1059842e-08,  0.2770890e-09; ...
             0.2117582e-21,   0.6708473e-26,  0.5652526e-10,  0.6541195e-13,  0.8709806e-14, -0.4430849e-22; ...
             0.5293956e-21,   0.1091152e-25,  0.6541195e-13,  0.5749008e-10,  0.7654984e-11, -0.7072071e-22; ...
             0.3465237e-01,   0.1059842e-08,  0.8709806e-14,  0.7654984e-11,  0.2232401e-03,  0.5830381e-04; ...
             0.9048777e-02,   0.2770890e-09, -0.4430849e-22, -0.7072071e-22,  0.5830381e-04,  0.1522727e-04];

           eqcovtrace = [ ...
             0.5703960e-13,   0.2562410e-13,  0.2218139e-11,  0.6530893e-14, -0.1258488e-17,  0.6777803e-17; ...
             0.2562410e-13,   0.6348393e-13, -0.1958979e-11,  0.7398763e-14, -0.4223513e-17,  0.2274645e-16; ...
             0.2218139e-11,  -0.1958979e-11,  0.3562517e-09,  0.2727333e-15,  0.2373188e-12, -0.1278121e-11; ...
             0.6530893e-14,   0.7398763e-14,  0.2727333e-15,  0.1256922e-14,  0.4930381e-31, -0.1972152e-30; ...
            -0.1258488e-17,  -0.4223513e-17,  0.2373188e-12,  0.4930381e-31,  0.2292326e-13, -0.2061165e-16; ...
             0.6777803e-17,   0.2274645e-16, -0.1278121e-11, -0.1972152e-30, -0.2061165e-16,  0.2288388e-13];

           printcov( cartcov,'ct','m',anom );
           fprintf(1,'Check the xxxx (mat*transpose) of cartcov \n' );
           fprintf(1,'%16e%16e%16e%16e%16e%16e\n',(cartcov*cartcov')' );
%            printcov( classcovtrace,'cl','t',anom );
%            classcovt = covunits( classcovtrace,anom,'cl','m' );
%            printcov( classcovt,'cl','m',anom );
%            printcov( eqcovtrace,'eq','t',anom );
%            eqcovt = covunits( eqcovtrace,anom,'eq','m' );
%            printcov( eqcovt,'eq','m',anom );
        end

        if testnum == 4
           reci = [4364.51524926;    4748.17602940;    2430.20427647];
           veci = [   5.87962414;      -4.10294944;      -2.53527819];

           year = 2000;
           mon  =  12;
           day  =  14;
           hr   =   5;
           min  =  25;
           sec  =   3.461;
           dut1 =  0.10597;
           dat  = 32;
           xp   =  0.0;
           yp   =  0.0;
           lod  =  0.0;
           terms = 2;
           timezone= 0;
           order = 106;

           % ----- in from usstratcom is lower diagonal!!!
           eqcov = [ ...
             4.68914e-11,  1.60090e-11,  1.64731e-10, -4.38141e-16, -1.41195e-10,  1.09999e-11, -8.49933e-12; ...
             1.60090e-11,  1.06881e-11,  6.39732e-11,  6.53124e-16, -5.60554e-11,  2.56099e-11, -4.44240e-13; ...
             1.64731e-10,  6.39732e-11,  6.71336e-10,  1.04455e-14, -6.37507e-10,  9.40554e-11,  7.34847e-12; ...
            -4.38141e-16,  6.53124e-16,  1.04455e-14,  2.10897e-17,  1.59272e-14,  1.03081e-14,  1.80744e-13; ...
            -1.41195e-10, -5.60554e-11, -6.37507e-10,  1.59272e-14,  7.59844e-10, -1.05437e-10,  1.77857e-10; ...
             1.09999e-11,  2.56099e-11,  9.40554e-11,  1.03081e-14, -1.05437e-10,  1.86472e-10,  3.99631e-11; ...
            -8.49933e-12, -4.44240e-13,  7.34847e-12,  1.80744e-13,  1.77857e-10,  3.99631e-11,  4.67207e-06;];

            printcov( eqcovtrace,'eq','t',anom );
            eqcovt = covunits( eqcovtrace,anom,'eq','m' );
            printcov( eqcovt,'eq','m',anom );
        end % if testnum

        if testnum == 5
           reci = [10127.26750234;    6972.89492052;    4902.05501566];
           veci = [   -3.38546947;       0.81341143;      -5.70012872];

           year = 2000;
           mon  =  12;
           day  =   1;
           hr   =   4;
           min  =  39;
           sec  =   8.735;
           dut1 =  0.11974;
           dat  = 32;
           xp   =  0.0;
           yp   =  0.0;
           lod  =  0.0;
           terms = 2;
           timezone= 0;
           order = 106;

           eqcovtrace = [ ...
             1.79533e-09,  9.35041e-10, -1.99413e-09, -1.63427e-13,  1.05541e-10, -5.77137e-10,  5.44105e-10; ...
             9.35041e-10,  1.04209e-09, -5.91847e-10, -9.78886e-14, -4.02330e-11, -4.39526e-10,  6.10449e-10; ...
            -1.99413e-09, -5.91847e-10,  3.63062e-09,  2.62623e-13, -3.98481e-10,  9.56574e-10,  2.86557e-09; ...
            -1.63427e-13, -9.78886e-14,  2.62623e-13,  5.56451e-16,  5.18072e-14, -2.60304e-14, -2.27216e-12; ...
             1.05541e-10, -4.02330e-11, -3.98481e-10,  5.18072e-14,  1.00164e-09,  2.40555e-10, -1.10189e-10; ...
            -5.77137e-10, -4.39526e-10,  9.56574e-10, -2.60304e-14,  2.40555e-10,  2.05855e-09, -1.97222e-11; ...
             5.44105e-10,  6.10449e-10,  2.86557e-09, -2.27216e-12, -1.10189e-10, -1.97222e-11,  1.13395e-05;];

            printcov( eqcovtrace,'eq','t',anom );
            eqcovt = covunits( eqcovtrace,anom,'eq','m' );
            printcov( eqcovt,'eq','m',anom );
        end % if testnum

        
%         a = 6860.7631;
%         ecc = 0.0010640;
%         p = a*(1.0 - ecc^2);
%         incl = 90.0/rad; %97.65184/rad;
%         omega = 79.54701/rad;
%         argp = 83.86041/rad;
%         nu = 65.21303/rad;
%         m = 65.10238/rad;
%         arglat = 0.0;
%         truelon = 0.0;
%         lonper = 0.0;
%         [reci,veci] = coe2rv ( p,ecc,incl,omega,argp,nu,arglat,truelon,lonper );        

          
        fprintf(1,' ---------------------------- begin tests ------------------------- \n' );
        fprintf(1,' ---------------------------- begin tests ------------------------- ' );
        anomeq1 = 'mean';  % true, mean
        anomeq2 = 'a';     % a, n
        anomflt = 'latlon'; % latlon  radec

        [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
                = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

        % --- convert the eci state into the various other state formats (classical, equinoctial, etc) 
        [cartstate,classstate,flstate,eqstate, fr] = setcov(reci,veci, ...
                year,mon,day,hr,min,sec,dut1,dat,ttt,jdut1,lod,xp,yp,terms,'y',strcat(anomeq1,anomeq2),anomflt,ddpsi,ddeps);
            
        fprintf(1,'==================== do the sensitivity tests \n' );
        cartcov = [ ...
                    100.0,  1.0e-2, 1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
                    1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
                    1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4; ...
                    1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6; ...
                    1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6; ...
                    1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001];
              
        fprintf(1,'1.  Cartesian Covariance \n');
        printcov( cartcov,'ct','m',strcat(anomeq1,'a') );

        
        % test partials
%         % partial a wrt rx
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         magr = mag(reci);
%         delta = reci(1) * 0.00001;
%         reci(1) = reci(1) + delta;
%         magr1 = mag(reci);
%         [p,a1,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         p0 = (a-a1)/delta;
%         p1 = 2.0*a^2*reci(1) / magr^3;
%         fprintf(1,' a wrt rx %14.14f  %14.14f \n', p0, p1);
%         
%         % partial n wrt rx
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         magr = mag(reci);
%         n = sqrt(mu/a^3);
%         recin = reci;
%         delta = recin(1) * 0.00001;
%         recin(1) = recin(1) + delta;
%         magr1 = mag(recin);
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (recin,veci);
%         n1 = sqrt(mu/a^3);
%         p0 = (n-n1)/delta;
%         p2 = -3.0*n1*a*reci(1)/magr1^3; 
%         fprintf(1,' n wrt rx %14.14f  %14.14f \n', p0, p2);
%         
%         % partial a wrt vx
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         magr = mag(veci);
%         vecin = veci;
%         delta = vecin(1) * 0.00001;
%         vecin(1) = vecin(1) + delta;
%         magv1 = mag(vecin);
%         [p,a1,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,vecin);
%         p0 = (a-a1)/delta;
%         p1 = 2.0*veci(1) / (n^2*a);
%         fprintf(1,' a wrt vx %14.14f  %14.14f \n', p0, p1);
%         
%         % partial n wrt vx
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         magr = mag(veci);
%         n = sqrt(mu/a^3);
%         vecin = veci;
%         delta = vecin(1) * 0.00001;
%         vecin(1) = vecin(1) + delta;
%         magv1 = mag(vecin);
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,vecin);
%         n1 = sqrt(mu/a^3);
%         p0 = (n-n1)/delta;
%         p2 = -3.0*vecin(1)/(n*a^2); 
%         fprintf(1,' n wrt vx %14.14f  %14.14f \n', p0, p2);
%         
%         % partial n wrt vz
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         magr = mag(veci);
%         n = sqrt(mu/a^3);
%         vecin = veci;
%         delta = vecin(3) * 0.00001;
%         vecin(3) = vecin(3) + delta;
%         magv1 = mag(vecin);
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,vecin);
%         n1 = sqrt(mu/a^3);
%         p0 = (n-n1)/delta;
%         p2 = -3.0*vecin(3)/(n*a^2); 
%         fprintf(1,' n wrt vz %14.14f  %14.14f \n', p0, p2);
%         
%         % partial rx wrt a
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         delta = a * 0.00001;
%         a = a + delta;
%         p = a*(1-ecc^2);
%         [reci1, veci1] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);
%         p0 = (reci(1)-reci1(1))/delta;
%         p1 = reci(1) / a;
%         fprintf(1,' rx wrt a %14.14f  %14.14f \n', p0, p1);
%         % partial rx wrt n
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         n = sqrt(mu/a^3);
%         delta = n * 0.00001;
%         n = n + delta;
%         a = (mu/n^2)^(1/3);
%         p = a*(1-ecc^2);
%         [reci1, veci1] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);
%         p0 = (reci(1)-reci1(1))/delta;
%         p1 = -2*reci(1) / (3*n);
%         fprintf(1,' rx wrt n %14.14f  %14.14f \n', p0, p1);
% 
%         % partial vx wrt a
%          [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%          delta = a * 0.00001;
%          a = a + delta;
%          p = a*(1-ecc^2);
%          [reci1, veci1] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);
%          p0 = (veci(1)-veci1(1))/delta;
%          p1 = -veci(1) / (2*a);
%          fprintf(1,' vx wrt a %14.14f  %14.14f \n', p0, p1);
%          % partial vx wrt n
%          [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%          n = sqrt(mu/a^3);
%          delta = n * 0.00001;
%          n = n + delta;
%          a = (mu/n^2)^(1/3);
%          p = a*(1-ecc^2);
%          [reci1, veci1] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);
%          p0 = (veci(1)-veci1(1))/delta;
%          p1 = veci(1) / (3*n);
%          fprintf(1,' vx wrt n %14.14f  %14.14f \n', p0, p1);
% 
%         % partial vz wrt a
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         delta = a * 0.00001;
%         a1 = a + delta;
%         p = a1*(1-ecc^2);
%         [recin,vecin] = coe2rv ( p,ecc,incl,omega,argp,nu,0.0,0.0,0.0 );
%         n1 = sqrt(mu/a1^3);
%         p0 = (veci(3)-vecin(3))/delta;
%         p2 = -veci(3)/(2*a); 
%         fprintf(1,' vz wrt a %14.14f  %14.14f \n', p0, p2);
%         % partial vz wrt n
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         n = sqrt(mu/a^3);
%         delta = n * 0.00001;
%         n1 = n + delta;
%         a1 = (mu/n1^2)^(1/3);
%         p = a1*(1-ecc^2);
%         [recin,vecin] = coe2rv ( p,ecc,incl,omega,argp,nu,0.0,0.0,0.0 );
%         p0 = (veci(3)-vecin(3))/delta;
%         p2 = veci(3)/(3*n); 
%         fprintf(1,' vz wrt n %14.14f  %14.14f \n', p0, p2);
%         % partial lat wrt rx
%         %[p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
%         delta = flstate(2) * 0.0000001;
%         latgc1 = flstate(2) + delta;
%         magr = flstate(5);
%             recef(1) = magr*0.001*cos(latgc1)*cos(flstate(1));  % in km
%             recef(2) = magr*0.001*cos(latgc1)*sin(flstate(1));
%             recef(3) = magr*0.001*sin(latgc1);
%             vecef = [0; 0; 0];
%             aecef = [0;0;0];
%             [recin,vecin,a] = ecef2eci(recef',vecef,aecef,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);
%               
%         n = sqrt(mu/a^3);
%         a1 = (mu/n1^2)^(1/3);
%         p = a1*(1-ecc^2);
%         [recin,vecin] = coe2rv ( p,ecc,incl,omega,argp,nu,0.0,0.0,0.0 );
%         p0 = (reci(1)-recin(1))/delta;
%         p1 = -magr*sin(flstate(1))*cos(flstate(2));
%         p2 = -reci(2)/sqrt(reci(1)^2 + reci(2)^2); 
%         fprintf(1,' lat wrt rx %14.14f  %14.14f \n', p0, p2);
%  pause        

        % paper approach
        fprintf(1,'2.  Classical Covariance from Cartesian #1 above \n');
        [classco,tmct2cl] = covct2clnew( cartcov,cartstate,strcat(anomeq1,'a') );
        printcov( classco,'cl','m',strcat(anomeq1,'a') );

%    fprintf(1,'2a.  Classical Covariance from Cartesian #1 above OLD \n');
%    [classcoold,tm] = covct2cl( cartcov,cartstate,strcat(anomeq1,'a') );
%    printcov( classcoold,'cl','m',strcat(anomeq1,'a') );

   %strcat(anomeq1,anomeq2)

        fprintf(1,'3.  Equinoctial Covariance from Classical #2 above \n');
        [eqco, tmcl2eq]    = covcl2eq( classco,classstate,strcat(anomeq1,anomeq2), fr );
        printcov( eqco,'eq','m',strcat(anomeq1,anomeq2) );
        eqcov = eqco; % save for later
         
        fprintf(1,'4.  Classical Covariance from Equinoctial #3 above \n');
        [classco1,tmeq2cl] = coveq2cl( eqco,eqstate,strcat(anomeq1,anomeq2), fr );
        printcov( classco1,'cl','m',strcat(anomeq1,anomeq2) );
      
        fprintf(1,'5.  Equinoctial Covariance from Cartesian #1 above \n');
        [eqco, tmct2eq]   = covct2eq( cartcov,cartstate,strcat(anomeq1,anomeq2), fr );
        printcov( eqco,'eq','m',strcat(anomeq1,anomeq2) );

        fprintf(1,'6.  Cartesian Covariance from Equinoctial #5 above \n');
        % try with intermediate one
        [cartco2,tmeq2cti] = coveq2ct( eqco, eqstate, strcat(anomeq1,anomeq2), fr );
        printcov( cartco2,'ct','m',strcat(anomeq1,anomeq2) ); 

        fprintf(1,'7.  Cartesian Covariance from Equinoctial #3 above \n');
        % try with old saved one
        [cartco2,tmeq2ct] = coveq2ct( eqcov, eqstate, strcat(anomeq1,anomeq2), fr );
        printcov( cartco2,'ct','m',strcat(anomeq1,anomeq2) ); 
      
        fprintf(1,'8.  Cartesian Covariance from Classical #2 above \n');
        [cartcov1, tmcl2ct]   = covcl2ctnew( classco,classstate,strcat(anomeq1,'a') );
        printcov( cartcov1,'ct','m',strcat(anomeq1,'a') );

%    fprintf(1,'8a.  Cartesian Covariance from Classical #2 above OLD \n');
%    [cartcovold,tm] = covcl2ct( classco,classstate,strcat(anomeq1,'a') );
%    printcov( cartcovold,'ct','m',strcat(anomeq1,'a') );

%          if strcmp(anom1,'latlon') == 1  % 1 is true
%              aeci=[0;0;0];
%              [recef,vecef,aecef] = eci2ecef(reci,veci,aeci,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);
%              cartstate = [recef(1);recef(2);recef(3);vecef(1);vecef(2);vecef(3)];
%          end

        fprintf(1,'9.  Flight Covariance from Cartesian #1 above \n');
        [flcov, tmct2fl] = covct2fl( cartcov,cartstate, anomflt, ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps);
        if strcmp(anomflt,'latlon') == 1  % 1 is true
            printcov( flcov,'fl','m',strcat(anomeq1,anomeq2) );
        else
            printcov( flcov,'sp','m',strcat(anomeq1,anomeq2) );
        end

        fprintf(1,'10. Cartesian Covariance from Flight #9 above \n');
        [cartco, tmfl2ct] = covfl2ct( flcov,flstate,anomflt, ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps );
        printcov( cartco,'ct','m',strcat(anomeq1,anomeq2) );

        fprintf(1,'11. Classical Covariance from Flight #9 above \n');
        fprintf(1,'\n-------- tm fl2cl --------- \n');
        [classcova, tmfl2cl] = covfl2cltest( flcov,flstate,anomflt, ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps );
        printcov( classcova,'cl','m',strcat(anomeq1,anomeq2) );
        

        fprintf(1,'Check inverses of each tm  = inv(tm) * tm \n');
        fprintf(1,'-------- tm ct2cl ---------\n');
        printcov( tmct2cl,'tm','m',strcat(anomeq1,'a') );
        fprintf(1,'\n');
        printcov( inv(tmct2cl)*tmct2cl,'tm','m',strcat(anomeq1,'a') );
        fprintf(1,'-------- tm cl2eq ---------\n');
        printcov( tmcl2eq,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'\n');
        printcov( inv(tmcl2eq)*tmcl2eq,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'-------- tm eq2cl ---------\n');
        printcov( tmeq2cl,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'\n');
        fprintf(1,'-------- tm eq2cl from cl2eqinv ---------\n');
        printcov( inv(tmcl2eq),'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'\n');
%break;
        printcov( inv(tmeq2cl)*tmeq2cl,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'-------- tm ct2eq ---------\n');
        printcov( tmct2eq,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'\n');
        printcov( inv(tmct2eq)*tmct2eq,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'-------- tm eq2ct ---------\n');
        printcov( tmeq2ct,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'\n');
        printcov( inv(tmeq2ct)*tmeq2ct,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'-------- tm cl2ct ---------\n');
        printcov( tmcl2ct,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'\n');
        printcov( inv(tmcl2ct)*tmcl2ct,'tm','m',strcat(anomeq1,anomeq2) );
       

        fprintf(1,'Check consistency of both approaches tmct2cl-inv(tmcl2ct) diff pct over 1e-18 \n');
        fprintf(1,'-------- accuracy of tm comparing ct2cl and cl2ct --------- \n');
        tm1 = tmct2cl;
        tm2 = inv(tmcl2ct);
        for i=1:6
            for j = 1:6
                if (abs( tm1(i,j) - tm2(i,j) ) < small) || (abs(tm1(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tm1(i,j)-tm2(i,j)) / tm1(i,j));
                end;
            end;
        end;
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );
  
% temp testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        fprintf(1,'-------- tm cl2ct new ---------\n');
        printcov( tmcl2ct,'tm','m',strcat(anomeq1,'a') );
        fprintf(1,'-------- tm ct2cl new ---------\n');
        printcov( tmct2cl,'tm','m',strcat(anomeq1,'a') );

        fprintf(1,'2.temp testing!!  find original transformations \n');
        [classco, tmct2clO]  = covct2clnewO( cartcov,cartstate,strcat(anomeq1,'a') );
        [cartcov1, tmcl2ctO] = covcl2ctnewO( classco,classstate,strcat(anomeq1,'a') );

        fprintf(1,'-------- tm cl2ct old ---------\n');
        printcov( tmcl2ctO,'tm','m',strcat(anomeq1,'a') );
        fprintf(1,'-------- tm ct2cl old ---------\n');
        printcov( tmct2clO,'tm','m',strcat(anomeq1,'a') );
  
        fprintf(1,'Check consistency of both approaches tmct2cl-inv(tmcl2ct) diff pct over 1e-18 \n');
        fprintf(1,'-------- accuracy of tm comparing ct2cl and cl2ct --------- \n');
        tm1 = tmct2clO;
        tm2 = inv(tmcl2ctO);
        for i=1:6
            for j = 1:6
                if (abs( tm1(i,j) - tm2(i,j) ) < small) || (abs(tm1(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tm1(i,j)-tm2(i,j)) / tm1(i,j));
                end;
            end;
        end;
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );
% temp testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
pause;
        
        
        
        fprintf(1,'-------- accuracy of tm comparing cl2eq and eq2cl --------- \n');
        tm1 = tmcl2eq;
        tm2 = inv(tmeq2cl);
        for i=1:6
            for j = 1:6
                if (abs( tm1(i,j) - tm2(i,j) ) < small) || (abs(tm1(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tm1(i,j)-tm2(i,j)) / tm1(i,j));
                end;
            end;
        end;
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

        
        fprintf(1,'-------- accuracy of tm comparing ct2eq and eq2ct --------- \n');
        tm1 = tmct2eq;
        tm2 = inv(tmeq2ct);
        for i=1:6
            for j = 1:6
                if (abs( tm1(i,j) - tm2(i,j) ) < small) || (abs(tm1(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tm1(i,j)-tm2(i,j)) / tm1(i,j));
                end;
            end;
        end;
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );


        fprintf(1,'-------- accuracy of tm comparing ct2eq and eq2ctintermediate --------- \n');
        tm1 = tmct2eq;
        tm2 = inv(tmeq2cti);
        for i=1:6
            for j = 1:6
                if (abs( tm1(i,j) - tm2(i,j) ) < small) || (abs(tm1(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tm1(i,j)-tm2(i,j)) / tm1(i,j));
                end;
            end;
        end;
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

        
        fprintf(1,'-------- accuracy of tm comparing cl2fl and fl2cl --------- \n');
        tm1 = tmct2fl;
        tm2 = inv(tmfl2ct);
        for i=1:6
            for j = 1:6
                if (abs( tm1(i,j) - tm2(i,j) ) < small) || (abs(tm1(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tm1(i,j)-tm2(i,j)) / tm1(i,j));
                end;
            end;
        end;
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

        
        
        temp = tmcl2eq * tmct2cl;
        temp1 = tmcl2ct * tmeq2cl;
    
        fprintf(1,'\n-------- tm combined ct2cl, cl2eq --------- \n');
        printcov( temp,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'\n-------- tm ct2eq --------- \n');
        printcov( tmct2eq,'tm','m',strcat(anomeq1,anomeq2) );
        
        fprintf(1,'\n-------- tm combined eq2cl, cl2ct --------- \n');
        printcov( temp1,'tm','m',strcat(anomeq1,anomeq2) );
        fprintf(1,'\n-------- tm ct2eq --------- \n');
        printcov( tmeq2ct,'tm','m',strcat(anomeq1,anomeq2) );
        
        fprintf(1,'-------- accuracy of test tm ct2eq --------- \n');
        tm1 = temp;
        tm2 = tmct2eq;
        for i=1:6
            for j = 1:6
                if (abs( tm1(i,j) - tm2(i,j) ) < small) || (abs(tm1(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tm1(i,j)-tm2(i,j)) / tm1(i,j));
                end;
            end;
        end;
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

        fprintf(1,'-------- accuracy of test tm eq2ct --------- \n');
        tm1 = temp1;
        tm2 = tmeq2ct;
        for i=1:6
            for j = 1:6
                if (abs( tm1(i,j) - tm2(i,j) ) < small) || (abs(tm1(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tm1(i,j)-tm2(i,j)) / tm1(i,j));
                end;
            end;
        end;
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

     
pause;
        %anom = 'true';
        [classco1, tm] = covct2clnew( cartcov,cartstate,strcat(anomeq1,anomeq2) );
        printcov( classco1,'cl','m',strcat(anomeq1,anomeq2) );
        
        [cartcov1, tm]   = covcl2ct( classco1,classstate,strcat(anomeq1,anomeq2) );
        printcov( cartcov1,'ct','m',strcat(anomeq1,anomeq2) );

        
pause;
        
        fprintf(1,'===================== do the eq-ct-eq tests \n' );
        eqcov = [ ...
                    1.0e-14, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16; ...
                    1.0e-16, 1.0e-14, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16; ...
                    1.0e-16, 1.0e-16, 1.0e-14, 1.0e-16, 1.0e-16, 1.0e-16; ...
                    1.0e-16, 1.0e-16, 1.0e-16, 1.0e-19, 1.0e-16, 1.0e-16; ...
                    1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-14, 1.0e-16; ...
                    1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-14];
        printcov( eqcov,'eq','m',anom );
        [classco,tm] = coveq2cl( eqcov,eqstate,anom);
        [cartco]  = covcl2ct( classco,classstate,anom );
        [classco,tm] = covct2cl( cartco,cartstate,anom );
        [eqco]    = covcl2eq( classco,classstate,anom );
        printcov( eqco,'eq','m',anom );
        for i=1:6
            for j = 1:6
                if (abs( eqcov(i,j)-eqco(i,j) ) < small) || (abs(eqcov(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (eqcov(i,j)-eqco(i,j)) / eqcov(i,j));
                  end;
            end;
        end;
        fprintf(1,'-------- accuracy --------- \n');
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );


        % ---------------- reset original cartcov ----------------------------
        cartcov = [ ...
                    1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
                    0.0, 1.0, 0.0, 0.0, 0.0, 0.0; ...
                    0.0, 0.0, 1.0, 0.0, 0.0, 0.0; ...
                    0.0, 0.0, 0.0, 0.000001, 0.0, 0.0; ...
                    0.0, 0.0, 0.0, 0.0, 0.000001, 0.0; ...
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.000001];
        cartcov = [ ...
                    0.81, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8; ...
                    1.0e-8, 0.81, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8; ...
                    1.0e-8, 1.0e-8, 0.81, 1.0e-8, 1.0e-8, 1.0e-8; ...
                    1.0e-8, 1.0e-8, 1.0e-8, 0.000001, 1.0e-8, 1.0e-8; ...
                    1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 0.000001, 1.0e-8; ...
                    1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 0.000001];
        cartcov = [ ...
                    100.0,    1.0e-2, 1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
                    1.0e-2, 100.0,    1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4; ...
                    1.0e-2, 1.0e-2, 100.0,    1.0e-4,   1.0e-4,   1.0e-4; ...
                    1.0e-4, 1.0e-4, 1.0e-4, 0.0001, 1.0e-6,   1.0e-6; ...
                    1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001, 1.0e-6; ...
                    1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001];

        % -------------------------------------------------------------
        % ------------- start covariance conversion tests -------------
        fprintf(1,'anomaly is %s \n',anom );
        fprintf(1,'initial cartesian covariance \n');
        printcov( cartcov,'ct','m',anom );

        fprintf(1,' -------- classical covariance conversions --------- \n');
        fprintf(1,' -----  cartesian to classical covariance  --------- \n');
        [classcov,tm] = covct2cl(cartcov,cartstate,anom );
        printcov( classcov,'cl','m',anom );

        fprintf(1,'tm ct2cl \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',tm' );
        fprintf(1,'tm ct2cl inv \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',inv(tm)' );
        tmold = inv(tm);

%        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
%        fprintf(1,'%16.12f%16.12f%16.12f%16.12f%16.12f%16.12f\n',(tm * tm')' );

%        fprintf(1,'Check the xxxx(transpose*mat) of classcov \n' );
%        fprintf(1,'%16.12f%16.12f%16.12f%16.12f%16.12f%16.12f\n',orth(classcov)' );
%        fprintf(1,'Check the xxxx(mat*transpose) of cartcov \n' );
%        fprintf(1,'%16.12f%16.12f%16.12f%16.12f%16.12f%16.12f\n',orth(cartcov)' );

%        % --- now check out what the scaling should be beteween
%        for i=1:6
%          for j = 1:6
%            diffm(i,j) = classcovtrace(i,j)/classcov(i,j);
%           end;
%         end;
%        fprintf(1,'scaling factors - should be \n');
%        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',diffm' );
%        for i=1:6
%           for j = 1:6
%                if (abs( classcovt(i,j)-classcov(i,j) ) < small) | (abs(classcovt(i,j)) < small)
%                    diffm(i,j) = 0.0;
%                  else
%                    diffm(i,j) = 100.0*( (classcovt(i,j)-classcov(i,j)) / classcovt(i,j));
%                  end;
%              end;
%          end;
%
%        % diff just to trace output
%        fprintf(1,'-------- accuracy to trace --------- \n');
%        fprintf(1,'pct differences if over %4e \n', small);
%        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

        fprintf(1,'----- classical to cartesian covariance  ---------- \n');
        [cartco,tm] = covcl2ct( classcov,classstate,anom );
        printcov( cartco,'ct','m',anom );

        fprintf(1,'tm cl2ct \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',tm' );

%        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
%        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',(tm*tm')' );

        for i=1:6
            for j = 1:6
                if (abs( cartcov(i,j)-cartco(i,j) ) < small) || (abs(cartcov(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (cartcov(i,j)-cartco(i,j)) / cartcov(i,j));
                  end;
              end;
          end;
        fprintf(1,'-------- accuracy --------- \n');
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

        for i=1:6
            for j = 1:6
                if (abs( tmold(i,j)-tm(i,j) ) < small) || (abs(tmold(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tmold(i,j)-tm(i,j)) / tmold(i,j));
                  end;
              end;
          end;
        fprintf(1,'-------- tm accuracy --------- \n');
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );
pause;
        fprintf(1,' -------- flight covariance conversions --------- \n');
        fprintf(1,' -----  cartesian to flight covariance  --------- \n');
        [flcov,tm] = covct2fl( cartcov,cartstate);
        printcov( flcov,'fl','m',anom );

        fprintf(1,'tm ct2fl \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',tm' );
        fprintf(1,'tm ct2fl inv \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',inv(tm)' );
        tmold = inv(tm);

%        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
%        fprintf(1,'%16.12f%16.12f%16.12f%16.12f%16.12f%16.12f\n',(tm * tm')' );

%        fprintf(1,'Check the xxxx(transpose*mat) of flcov \n' );
%        fprintf(1,'%16.12f%16.12f%16.12f%16.12f%16.12f%16.12f\n',(flcov' * flcov)' );

        fprintf(1,'----- flight to cartesian covariance \n');
        [cartco,tm] = covfl2ct( flcov,flstate,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps );
        printcov( cartco,'ct','m',anom );

        fprintf(1,'tm  fl2ct \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',tm' );

%        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
%        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',(tm*tm')' );
        for i=1:6
            for j = 1:6
                if (abs( cartcov(i,j)-cartco(i,j) ) < small) | (abs(cartcov(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (cartcov(i,j)-cartco(i,j)) / cartcov(i,j));
                  end;
              end;
          end;
        fprintf(1,'-------- accuracy --------- \n');
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

        for i=1:6
            for j = 1:6
                if (abs( tmold(i,j)-tm(i,j) ) < small) | (abs(tmold(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tmold(i,j)-tm(i,j)) / tmold(i,j));
                  end;
              end;
          end;
        fprintf(1,'-------- tm accuracy --------- \n');
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );
pause;

        fprintf(1,' -------- equinoctial covariance conversions --------- \n');
        fprintf(1,' -----  classical to equinoctial covariance ---------- \n');
        [eqcov,tm] = covcl2eq ( classcov,classstate,anom);
        printcov( eqcov,'eq','m',anom );

        fprintf(1,'tm cl2eq \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',tm' );
        fprintf(1,'tm cl2eq inv \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',inv(tm)' );
        tmold = inv(tm);

%        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
%        fprintf(1,'%16.12f%16.12f%16.12f%16.12f%16.12f%16.12f\n',(tm * tm')' );

%        fprintf(1,'Check the xxxx(transpose*mat) of eqcov \n' );
%        fprintf(1,'%16.12f%16.12f%16.12f%16.12f%16.12f%16.12f\n',(eqcov' * eqcov)' );

        fprintf(1,'----- equinoctial to classical covariance \n');
        [classco,tm] = coveq2cl( eqcov,eqstate,anom);
        printcov( classco,'cl','m',anom );

        fprintf(1,'tm eq2cl \n' );
        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',tm' );

%        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
%        fprintf(1,'%16e%16e%16e%16e%16e%16e\n',(tm*tm')' );

        for i=1:6
            for j = 1:6
                if (abs( classcov(i,j)-classco(i,j) ) < small) | (abs(classcov(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (classcov(i,j)-classco(i,j)) / classcov(i,j));
                  end;
              end;
          end;
        fprintf(1,'-------- accuracy --------- \n');
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );

        for i=1:6
            for j = 1:6
                if (abs( tmold(i,j)-tm(i,j) ) < small) | (abs(tmold(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (tmold(i,j)-tm(i,j)) / tmold(i,j));
                  end;
              end;
          end;
        fprintf(1,'-------- accuracy tm --------- \n');
        fprintf(1,'pct differences if over %4e \n',small);
        fprintf(1,'%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n',diffm' );



pause;


        % ------------------------ check nav tests -----------------------------------
        small = 1.0e-18;
        doall = 'n';
%        fprintf(1,'year mon day hms  magr  magv  r sig  vsig  max diag  max mat \n');
        fprintf(1,'      ecc         incl        maxdiag        maxdiff        magr        r sig   \n');
        for testnumbb = 1:12  % total is 122
%        for testnumb = 1:122
            if testnumbb == 1
                testnumb = 3;
                satnum =   107;
              end;
            if testnumbb == 2
                testnumb = 16;
                satnum =    11;
              end;
            if testnumbb == 3
                testnumb = 28;
                satnum = 16609;
              end;
            if testnumbb == 4
                testnumb = 42;
                satnum = 20052;
              end;
            if testnumbb == 5
                testnumb = 52;
                satnum = 21867;
              end;
            if testnumbb == 6
                testnumb = 64;
                satnum = 23019;
              end;
            if testnumbb == 7
                testnumb = 76;
                satnum = 24780;
              end;
            if testnumbb == 8
                testnumb = 87;
                satnum = 25013;
              end;
            if testnumbb == 9
                testnumb = 97;
                satnum = 25544;
              end;
            if testnumbb == 10
                testnumb = 101;
                satnum = 25634;
              end;
            if testnumbb == 11
                testnumb = 110;
                satnum = 26354;
              end;
            if testnumbb == 12
                testnumb = 113;
                satnum = 26405;
              end;
            testcove;


        fprintf(1,'year %5i ',year);
        fprintf(1,' mon %4i ',mon);
        fprintf(1,' day %3i ',day);
        fprintf(1,'%3i:%2i:%8.6f\n ',hr,min,sec );
        fprintf(1,'dut1 %8.6f s',dut1);
        fprintf(1,' dat %3i s',dat);
        fprintf(1,' xp %8.6f "',xp);
        fprintf(1,' yp %8.6f "',yp);
        fprintf(1,' lod %8.6f s\n',lod);
        fprintf(1,' ddpsi %8.6f " ddeps  %8.6f\n',ddpsi,ddeps);
%        fprintf(1,'order %3i  eqeterms %31  opt %3s \n',order,eqeterms,opt );
        fprintf(1,'units are km and km/s and km/s2\n' );
        fprintf(1,'eci%14.7f%14.7f%14.7f',reci );
        fprintf(1,' v %14.9f%14.9f%14.9f\n',veci );

          [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
                  = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

        fprintf(1,'ut1 %8.6f tut1 %16.12f jdut1 %18.11f\n',ut1,tut1,jdut1 );
        fprintf(1,'utc %8.6f\n',utc );
        fprintf(1,'tai %8.6f\n',tai );
        fprintf(1,'tt  %8.6f ttt  %16.12f jdtt  %18.11f\n',tt,ttt,jdtt );
        fprintf(1,'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n',tdb,ttdb,jdtdb );

            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (reci,veci);
            fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
            fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n',...
                    p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
                    arglat*rad,truelon*rad,lonper*rad );


            [cartstate,classstate,flstate,eqstate] = setcov(reci,veci, ...
            year,mon,day,hr,min,sec,dut1,dat,ttt,jdut1,lod,xp,yp,terms,'n',anom,ddpsi,ddeps);




            fprintf(1,'new case --------------\n');

           for i=1:6
               eqcov(i,4) = eqcov(i,4)/tusec;
             end;
           for i=1:6
               eqcov(4,i) = eqcov(4,i)/tusec;
             end;

            % --------------------- write out input data --------------------------
            if doall ~= 'y'
%                printcov( eqcov,'eq','m',anom );
              end;
            printcov( eqcov,'eq','m',anom );

            % ------ do conversions
            [classcov,tm] = coveq2cl( eqcov,eqstate,anom);
            printcov( classcov,'cl','m',anom );

            [cartcov,tm] = covcl2ct( classcov,classstate,anom );
            printcov( cartcov,'ct','m',anom );

            [covoprsw,tm] = covct2o1(cartcov,cartstate);
            fprintf(1,'rsw\n');
            printcov( covoprsw,'ct','m',anom );
            temt = covoprsw;
            fprintf(1,'rsw %11.3f  %11.3f  %11.3f %11.3f \n',sqrt(temt(1,1)), ...
                       sqrt(temt(2,2)),sqrt(temt(3,3)), ...
                       sqrt(temt(1,1) + temt(2,2) + temt(3,3)) );

            [covopntw,tm] = covct2o2(cartcov,cartstate);
            fprintf(1,'ntw\n');
            printcov( covopntw,'ct','m',anom );
            temt = covopntw;
            fprintf(1,'rsw %11.3f  %11.3f  %11.3f %11.3f \n',sqrt(temt(1,1)), ...
                       sqrt(temt(2,2)),sqrt(temt(3,3)), ...
                       sqrt(temt(1,1) + temt(2,2) + temt(3,3)) );

            [classco,tm] = covct2cl(cartcov,cartstate,anom );
            printcov( classco,'cl','m',anom );

            [eqco,tm] = covcl2eq ( classco,classstate,anom);
            if doall ~= 'y'
%                printcov( eqco,'eq','m',anom );
              end;
            printcov( eqco,'eq','m',anom );

            for i=1:6
              for j = 1:6
                if (abs( eqcov(i,j)-eqco(i,j) ) < small) | (abs(eqcov(i,j)) < small)
                    diffm(i,j) = 0.0;
                  else
                    diffm(i,j) = 100.0*( (eqcov(i,j)-eqco(i,j)) / eqcov(i,j));
                  end;
               end;
             end;
            if doall ~= 'y'
%                fprintf(1,'accuracy \n');
%                fprintf(1,'%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f\n',diffm' );
              end;
%            fprintf(1,'accuracy \n');
%            fprintf(1,'%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f\n',diffm' );

            magr = sqrt( cartcov(1,1) + cartcov(2,2) + cartcov(3,3) );
            magv = sqrt( cartcov(4,4) + cartcov(5,5) + cartcov(6,6) );

            magrs = sqrt( possigma(1)^2 + possigma(2)^2 + possigma(3)^2 )*1000.0;
            magvs = sqrt( velsigma(1)^2 + velsigma(2)^2 + velsigma(3)^2 )*1000.0;

            if doall ~= 'y'
                fprintf(1,'position mag %11.5f velocity %11.5f \n',magr,magv );
                fprintf(1,'position sig from msg %11.5f velocity %11.5f \n',magrs,magvs );
                fprintf(1,'-------------------------------------------------\n' );
              end;

            md = max(diag(abs(diffm)));
            mm = max(max(abs(diffm)));

            if doall ~= 'y'
%                fprintf(1,'%3i %5i %3i %3i %3i:%2i:%8.6f %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n', ...
%                          testnumb,year,mon,day,hr,min,sec,magr,magv,magrs,magvs,md,mm );
%                if ( md > 0.1 ) | ( mm > 0.1 )
%                    fprintf(1,'input equinoctial covariance \n');
%                    printcov( eqcov,'eq','m',anom );
%                    fprintf(1,'calculated equinoctial covariance \n');
%                    printcov( eqco,'eq','m',anom );
%                  end
              end;
            fprintf(1,'%11.6f   %11.6f   %14.6f %14.6f %14.6f %11.6f %3i \n', ...
                    classstate(2),classstate(3)*rad,md,mm,magr,magrs,testnumb );

          end;



