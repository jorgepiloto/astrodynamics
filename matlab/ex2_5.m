%     -----------------------------------------------------------------
%
%                              Ex2_5.m
%
%  this file demonstrates example 2-5. it also includes some sttressing
%  cases forthe coe and rv conversions for all orbit types. 
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

            rad = 180.0 /pi;

            fprintf(1,'coe test ----------------------------\n' );
            r=[ 6524.834;6862.875;6448.296];
            v=[ 4.901327;5.533756;-1.976341];

            fprintf(1,'start %15.9f %15.9f %15.9f',r );
            fprintf(1,' v %15.10f %15.10f %15.10f\n',v );
            % --------  coe2rv       - classical elements to posisiotn and velocity
            % --------  rv2coe       - position and velocity vectors to classical elements
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeS (r,v);
            fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
            fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
                    p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
                    arglat*rad,truelon*rad,lonper*rad );
pause;
% test various combinations of coe and rv
        fprintf(1,'coe tests ----------------------------\n' );
        for i = 1:20
            if i == 1
                r = r;
                v = v;
              end
            if i == 2
                fprintf(1,'coe test ----------------------------\n' );
                r=[ 6524.834;6862.875;6448.296];
                v=[ 4.901327;5.533756;-1.976341];
              end

            % ------- elliptical orbit tests -------------------
            if i == 3
                fprintf(1,'coe test elliptical ----------------------------\n' );
                r=[ 1.1372844; -1.0534274; -0.8550194]*6378.137;
                v=[0.6510489;  0.4521008;  0.0381088]*7.905366149846;
              end
            if i == 4
                fprintf(1,'coe test elliptical ----------------------------\n' );
                r=[  1.0561942;-0.8950922;-0.0823703]*6378.137;;
                v=[  -0.5981066;-0.6293575; 0.1468194]*7.905366149846;
              end

            % ------- circular inclined orbit tests -------------------
            if i == 5
                fprintf(1,'coe test near circular inclined ----------------------------\n' );
                r=[ -0.4222777; 1.0078857; 0.7041832]*6378.137;
                v=[  -0.5002738;-0.5415267; 0.4750788]*7.905366149846;
              end
            if i == 6
                fprintf(1,'coe test near circular inclined ----------------------------\n' );
                r=[ -0.7309361;-0.6794646;-0.8331183]*6378.137;
                v=[  -0.6724131; 0.0341802; 0.5620652]*7.905366149846;
              end

            if i == 7 % -- CI u = 45 deg
                fprintf(1,'coe test circular inclined ----------------------------\n' );
                r = [-2693.34555010128  6428.43425355863  4491.37782050409];
                v = [   -3.95484712246016  -4.28096585381370  3.75567104538731];
              end
            if i == 8 % -- CI u = 315 deg
                fprintf(1,'coe test circular inclined ----------------------------\n' );
                r = [-7079.68834483379;  3167.87718823353; -2931.53867301568];
                v = [    1.77608080328182;  6.23770933190509; 2.45134017949138];
              end

            % ------- elliptical equatorial orbit tests -------------------
            if i == 9
                fprintf(1,'coe test elliptical near equatorial ----------------------------\n' );
                r=[ 21648.6109280739; -14058.7723188698; -0.0003598029];
                v=[ 2.16378060719980; 3.32694348486311; 0.00000004164788 ];
              end
            if i == 10
                fprintf(1,'coe test elliptical near equatorial ----------------------------\n' );
                r=[  7546.9914487222;  24685.1032834356; -0.0003598029];
                v=[ 3.79607016047138; -1.15773520476223; 0.00000004164788 ];
              end

            if i == 11 % -- EE w = 20 deg
                fprintf(1,'coe test elliptical equatorial ----------------------------\n' );
                r = [-22739.1086596208  -22739.1086596208     0.0];
                v = [    2.48514004188565  -2.02004112073465  0.0];
              end
            if i == 12 % -- EE w = 240 deg
                fprintf(1,'coe test elliptical equatorial ----------------------------\n' );
                r = [ 28242.3662822040    2470.8868808397    0.0];
                v = [    0.66575215057746  -3.62533022188304  0.0];
              end

            % ------- circular equatorial orbit tests -------------------
            if i == 13
                fprintf(1,'coe test circular near equatorial ----------------------------\n' );
                r=[ -2547.3697454933; 14446.8517254604; 0.000 ];
                v=[  -5.13345156333487; -0.90516601477599; 0.00000090977789 ];
              end
            if i == 14
                fprintf(1,'coe test circular near equatorial ----------------------------\n' );
                r=[  7334.858850000; -12704.3481945462;   0.000 ];
                v=[  -4.51428154312046; -2.60632166411836; 0.00000090977789 ];
              end

            if i == 15 % -- CE l = 65 deg
                fprintf(1,'coe test circular equatorial ----------------------------\n' );
                r = [ 6199.6905946008; 13295.2793851394;      0.0];
                v = [ -4.72425923942564; 2.20295826245369;    0.0];
              end
            if i == 16 % -- CE l = 65 deg i = 180 deg
                fprintf(1,'coe test circular equatorial ----------------------------\n' );
                r = [ 6199.6905946008; -13295.2793851394;      0.0];
                v = [ -4.72425923942564; -2.20295826245369;    0.0];
              end

            % ------- parabolic orbit tests -------------------
            if i == 17
                fprintf(1,'coe test parabolic ----------------------------\n' );
                r=[  0.5916109;-1.2889359;-0.3738343]*6378.137;
                v=[   1.1486347;-0.0808249;-0.1942733]*7.905366149846;
              end

            if i == 18
                fprintf(1,'coe test parabolic ----------------------------\n' );
                r=[-1.0343646; -0.4814891;  0.1735524]*6378.137;
                v=[ 0.1322278; 0.7785322; 1.0532856  ]*7.905366149846;
              end

            if i == 19
                fprintf(1,'coe test hyperbolic ---------------------------\n' );
                r=[0.9163903; 0.7005747; -1.3909623  ]*6378.137;
                v=[0.1712704; 1.1036199; -0.3810377  ]*7.905366149846;
              end

            if i == 20
                fprintf(1,'coe test hyperbolic ---------------------------\n' );
                r=[12.3160223; -7.0604653; -3.7883759]*6378.137;
                v=[-0.5902725; 0.2165037; 0.1628339  ]*7.905366149846;
              end

            fprintf(1,'start %15.9f %15.9f %15.9f',r );
            fprintf(1,' v  %15.10f %15.10f %15.10f\n',v );
            % --------  coe2rv       - classical elements to posisiotn and velocity
            % --------  rv2coe       - position and velocity vectors to classical elements
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r,v);
            fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
            fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
                    p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
                    arglat*rad,truelon*rad,lonper*rad );
            [r,v] = coe2rv(p,ecc,incl,omega,argp,nu,arglat,truelon,lonper);
            fprintf(1,'r     %15.9f %15.9f %15.9f',r );
            fprintf(1,' v  %15.10f %15.10f %15.10f\n',v );

        end; % for

