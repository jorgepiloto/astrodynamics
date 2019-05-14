    %
    % test the hills conversion programs
    % this does the intiial checkouts
    %  use testhillSTK for te longer runs with perturbed values, etc
    %
    % Uses:
    %  EQCM_to_ECI_RTN_sal
    %      f_ECI_to_RTN_sal
    %      newtonnu, newtone
    %      inverselliptic2
    %      elliptic12
    %  ECI_to_EQCM_RTN_sal
    %
    %
    %
    %
    % STK needs to be running (testhillsc)
    %
    %
    %

    fprintf(1,'------------------------------------------------ initial accuracy ----------------------------- \n');
    % now test the ability to convert eci - hills and back
    constmath;
    constastro;

    casenumo = 1;  % satellite test to conduct. paper used case 1 and case 8
    casetest = 7;  % 100 m in each axis, .01 m/s each axis
    ang_step = 0.00000001;

    % case where you read in the int and tgt ephemerides
    ropt = 'm'; % matlab or stk

    outfilehill = fopen(strcat('d:/STKFiles Educational Files/Hills/thillc', int2str(casenumo),'n', int2str(casetest), '.out'), 'wt');      % rewrite each time

    % -----------------------------------------------------------------------------------------------
    % -------------------------------- do initial accuracy checks fwd, back, ------------------------
    % -----------------------------------------------------------------------------------------------
    % --- this one also sets the case for the rest of the program here...
    for ktr = 1:10  % for closure checks
        % set last one to be the case of interest
        if ktr < 10
            casenum = ktr;
        else
            casenum = casenumo;
        end
        switch casenum
            case 1  % baseline circular leo
                rtgteci(1)= (6378.137 + 500.0); % km circular example
                rtgteci(2)= 0.0;
                rtgteci(3)= 0.0;
                magrt = mag( rtgteci );
                vtgteci(1)= 0.0;
                vtgteci(2)= sqrt(mu/magrt);   % make perfect circular orbit
                vtgteci(3)= 0.0;
                circ = 'y';
                rtgteci = rtgteci';  % get in col format
                vtgteci = vtgteci';

                %  previous doesn't work in STK92 - fixed in 9.3 for smaller inclination
                a = 6378.137 + 500.0;
                ecc = 0.0;
                p = a * (1.0 - ecc*ecc);
                [rtgteci,vtgteci] = coe2rv(p,0.0,0.001/rad,0.0,0.0,0.0, 0.0,0.0,0.0);
                circ = 'y';
                fprintf(1,'rtgt = [%20.13f %20.13f %20.13f]; \n vtgt = [%20.13f %20.13f %20.13f]; \n', ...
                    rtgteci(1),rtgteci(2),rtgteci(3),vtgteci(1),vtgteci(2),vtgteci(3) );
            case 2  % two-body gps equatorial
                a = 26500.0;
                ecc = 0.0;
                p = a * (1.0 - ecc*ecc);
                [rtgteci,vtgteci] = coe2rv(p,0.0,0.001/rad,0.0,0.0,0.0, 0.0,0.0,0.0);
                circ = 'y';
            case 3  % two-body gps inclined
                a = 26500.0;
                ecc = 0.0;
                p = a * (1.0 - ecc*ecc);
                [rtgteci,vtgteci] = coe2rv(p,0.0,55.0/rad,0.0,0.0,0.0, 0.0,0.0,0.0);
                circ = 'y';
            case 4  % two-body geo
                a = 42164.0;
                ecc = 0.0;
                p = a * (1.0 - ecc*ecc);
                [rtgteci,vtgteci] = coe2rv(p,0.0,0.001/rad,0.0,0.0,0.0, 0.0,0.0,0.0);
                circ = 'y';
            case 5  % LEO some eccentricity
                rtgteci = [-605.79043080 -5870.23040700 3493.05200400]';   % some eccenticity and inclination
                vtgteci = [-1.568251615 -3.702348353 -6.479484915]';
                circ = 'n';
            case 6  % gps inclined, slight ecc
                rtgteci = [-761.24075519  22699.40899449  13644.17603200]';
                vtgteci = [-2.630715586  1.396719096  -2.491891000]';
                circ = 'n';
            case 7  % geo equaotrial ecc
                rtgteci = [-40588.150362 -11462.167028 27.147649  ]';
                vtgteci = [0.834787457 -2.958305691 -0.001173016  ]';
                circ = 'n';
            case 8  % LEO higher, inclined, ecc
                rtgteci = [-5551.89864600	-2563.04969600	3257.75616500]';  % satellite 11
                vtgteci = [2.149073000	-7.539457000	-2.185709000]';
                circ = 'n';
            case 9  % heo high ecc
                rtgteci = [9668.14551571	6999.07240705	4041.43303674]';
                vtgteci = [-3.652923060	0.649665190	-5.821235160]';
                circ = 'n';
            case 10  % LEO higher, inclined, ecc
                a = 8164.7188;
                ecc = 0.08;
                p = a * (1.0 - ecc*ecc);
                [rtgteci,vtgteci] = coe2rv(p,ecc,32.861/rad,0.0, 0.0,0.0, 0.0,0.0,0.0);
                circ = 'n';
        end; % switch

        x  = 10.0;  % in m
        y  = 10.0;  %100.0;  % in m
        z  = 10.0;  %100.0;  % in m
        xd = 0.0100;  %0.10;
        yd = 0.0100;  %0.10;
        zd = 0.0100;  %0.10;

        hro1(1)=  x/1000.0; % km
        hro1(2)=  y/1000.0;
        hro1(3)=  z/1000.0;
        hrokm  = hro1'; % get in col format
        hvo1(1)= xd/1000.0; % km/s
        hvo1(2)= yd/1000.0;
        hvo1(3)= zd/1000.0;
        hvokm  = hvo1'; % get in col format

        [rintecix, vintecix]   = EQCM_to_ECI_RTN_sal( rtgteci, vtgteci, hrokm, hvokm );     % find int pos eqcm
        [rhillx, vhillx]      = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci, rintecix,vintecix, outfilehill ); % find how much this would be for true hills

        dr = rhillx - hrokm;
        dv = vhillx - hvokm;

        mdr = mag(dr*1000000);  % into mm
        mdv = mag(dv*1000000);
        fprintf(1,'dr m         %20.13f %20.13f %20.13f %15.8f mm %20.13f %20.13f %20.13f %15.8f mm/s \n', ...
            dr(1)*1000,dr(2)*1000,dr(3)*1000,mdr,dv(1)*1000,dv(2)*1000,dv(3)*1000,mdv );

        [p,a,ecc,incl,omega,argpx1,nux1,m,arglat,truelon,lonper ] = rv2coe (rtgteci, vtgteci);
        fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
            p,a,ecc,incl*rad,omega*rad,argpx1*rad,nux1*rad,m*rad, ...
            arglat*rad,truelon*rad,lonper*rad );

        [p,a,ecc,incl,omega,argpx1,nux1,m,arglat,truelon,lonper ] = rv2coe (rintecix, vintecix);
        fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
            p,a,ecc,incl*rad,omega*rad,argpx1*rad,nux1*rad,m*rad, ...
            arglat*rad,truelon*rad,lonper*rad );

    end;  % for ktr through 1-10 validation
    pause;

    % -----------------------------------------------------------------------------------------------
    % ------------------- check various positions to determine if the veocity is correct ------------
    % -----------------------------------------------------------------------------------------------
    for ktr = 1:6  % for initial setup
        if ktr == 1
            x = 100.0;  % in m
        else
            x = 0.0;
        end
        if ktr == 2
            y = 100.0;
        else
            y = 0.0;
        end
        if ktr == 3
            z = 100.0;
        else
            z = 0.0;
        end
        if ktr == 4
            xd = 0.01;
        else
            xd = 0.0;
        end
        if ktr == 5
            yd = 0.01;
        else
            yd = 0.0;
        end
        if ktr == 6
            zd = 0.01;
        else
            zd = 0.0;
        end

        hro1(1)=  x/1000.0; % km
        hro1(2)=  y/1000.0;
        hro1(3)=  z/1000.0;
        hrokm  = hro1'; % get in col format
        hvo1(1)= xd/1000.0; % km/s
        hvo1(2)= yd/1000.0;
        hvo1(3)= zd/1000.0;
        hvokm  = hvo1'; % get in col format

        [rintecisal,vintecisal]   = EQCM_to_ECI_RTN_sal( rtgteci,vtgteci,hrokm,hvokm );     % find int pos eqcm

        fprintf(1,'hillsin       %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n', ...
            hrokm(1)*1000,hrokm(2)*1000,hrokm(3)*1000,hvokm(1)*1000,hvokm(2)*1000,hvokm(3)*1000 );

        mmagrti = mag(rtgteci);
        mmagvti = mag(vtgteci);
        fprintf(1,'rtgteci        %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f %15.8f %15.8f \n', ...
            rtgteci(1),rtgteci(2),rtgteci(3),vtgteci(1),vtgteci(2),vtgteci(3),mmagrti,mmagvti );

        mmagri = mag(rintecisal);
        mmagvi = mag(vintecisal);
        fprintf(1,'rinteci        %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f  %15.8f %15.8f \n', ...
            rintecisal(1),rintecisal(2),rintecisal(3),vintecisal(1),vintecisal(2),vintecisal(3),mmagri, mmagvi );
    end  % for ktr 1-6


    pause;

    % -----------------------------------------------------------------------------------------------
    % -------------------------- perturb each one, or all of the vector components ------------------
    % -----------------------------------------------------------------------------------------------
    fprintf(1,'\n\n------------------------------- perturb each one ---------------------------------------- \n');
    % initialize
    x  = 0.0;  % in m
    y  = 0.0;  %100.0;  % in m
    z  = 0.0;  %100.0;  % in m
    xd = 0.000;  %0.10;
    yd = 0.000;  %0.10;
    zd = 0.000;  %0.10;
    fid = 2;
    for ktr = 1:6  % run these 6 cases...
        if ktr == 1
            x = 100.0;  % in m
        else
            x = 0.0;
        end
        if ktr == 2
            y = 100.0;
        else
            y = 0.0;
        end
        if ktr == 3
            z = 100.0;
        else
            z = 0.0;
        end
        if ktr == 4
            xd = 1;
        else
            xd = 0.0;
        end
        if ktr == 5
            yd = 1;
        else
            yd = 0.0;
        end
        if ktr == 6
            zd = 1;
        else
            zd = 0.0;
        end

        hro1(1)=  x/1000.0; % km
        hro1(2)=  y/1000.0;
        hro1(3)=  z/1000.0;
        hrokm = hro1'; % get in col format
        hvo1(1)= xd/1000.0; % km/s
        hvo1(2)= yd/1000.0;
        hvo1(3)= zd/1000.0;
        hvokm = hvo1'; % get in col format
        % reset this!! so m for hills call!!
        hro(1)=  x; % m
        hro(2)=  y;
        hro(3)=  z;
        hvo(1)= xd; % m/s
        hvo(2)= yd;
        hvo(3)= zd;

        [rintecisal,vintecisal]       = EQCM_to_ECI_RTN_sal( rtgteci,vtgteci,hrokm,hvokm );     % find int pos eqcm
        [rhill2, vhill2]              = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci,rintecisal,vintecisal, outfilehill ); % find how much this would be for true hills
        [rintecisal1,vintecisal1]     = EQCM_to_ECI_RTN_sal( rtgteci,vtgteci,rhill2,vhill2 );     % find int pos eqcm
        [rhillsal, vhillsal]          = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci,rintecisal1,vintecisal1, outfilehill ); % find how much this would be for true hills
        [rhillsal1, vhillsal1]        = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci,rintecisal,vintecisal, outfilehill ); % find how much this would be for true hills

        % ------------ set the vectors for future use in program
        rinteci(1) = rintecisal(1);
        rinteci(2) = rintecisal(2);
        rinteci(3) = rintecisal(3);
        vinteci(1) = vintecisal(1);
        vinteci(2) = vintecisal(2);
        vinteci(3) = vintecisal(3);

        % set initial values to compare and make sure STK doesn't have an error
        rintecio = rinteci;
        vintecio = vinteci;
        rtgtecio = rtgteci;
        vtgtecio = vtgteci;

        % --------------------- write out various transformation methods --- just at ktr = 8
        if (ktr <= 8)  %  || (ktr ==10)
            fprintf(1,'\nhills in ell   %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n', ...
                hro(1),hro(2),hro(3),hvo(1),hvo(2),hvo(3) );
            fprintf(1,'\nhills in       %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n', ...
                hrokm(1)*1000,hrokm(2)*1000,hrokm(3)*1000,hvokm(1)*1000,hvokm(2)*1000,hvokm(3)*1000 );

            mmagrti = mag(rtgteci);
            mmagvti = mag(vtgteci);
            fprintf(1,'rtgteci        %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n', ...
                rtgteci(1),rtgteci(2),rtgteci(3),mmagrti,vtgteci(1),vtgteci(2),vtgteci(3),mmagvti );

            fprintf(1,'rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n', ...
                rhill2(1)*1000,rhill2(2)*1000,rhill2(3)*1000,vhill2(1)*1000,vhill2(2)*1000,vhill2(3)*1000 );

            mmagrti = mag(rinteci);
            mmagvti = mag(vinteci);
            fprintf(1,'rinteci        %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n', ...
                rinteci(1),rinteci(2),rinteci(3),mmagrti,vinteci(1),vinteci(2),vinteci(3),mmagvti );

            mmagrti = mag(rintecisal);
            mmagvti = mag(vintecisal);
            fprintf(1,'rintecisal hil %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n', ...
                rintecisal(1),rintecisal(2),rintecisal(3),mmagrti,vintecisal(1),vintecisal(2),vintecisal(3),mmagvti );

            mmagrti = mag(rintecisal1);
            mmagvti = mag(vintecisal1);
            fprintf(1,'rintecisal1 eq %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n', ...
                rintecisal1(1),rintecisal1(2),rintecisal1(3),mmagrti,vintecisal1(1),vintecisal1(2),vintecisal1(3),mmagvti );

            fprintf(1,'hillssal hill  %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n', ...
                rhillsal(1)*1000,rhillsal(2)*1000,rhillsal(3)*1000,vhillsal(1)*1000,vhillsal(2)*1000,vhillsal(3)*1000 );
            fprintf(1,'hillssal1 eq   %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n', ...
                rhillsal1(1)*1000,rhillsal1(2)*1000,rhillsal1(3)*1000,vhillsal1(1)*1000,vhillsal1(2)*1000,vhillsal1(3)*1000 );

            fprintf(1,'==== sals codes ====\n');
            %---            [rintecisal,vintecisal]   = EQCM_to_ECI_NTW_sal( rtgteci*1000,vtgteci*1000,hro',hvo' );      % do in m~!!!find int pos eqcm
            %---            [rhill2, vhill2]          = ECI_to_EQCM_NTW_sal( rtgteci*1000,vtgteci*1000,rintecisal,vintecisal ); % find how much this would be for true hills
            [rintecisal, vintecisal]   = EQCM_to_ECI_RTN_sal( rtgteci,vtgteci,hro',hvo' );     % find int pos eqcm
            [rhill2, vhill2]      = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci,rintecisal,vintecisal, outfilehill ); % find how much this would be for true hills

            mmagrti = mag(rintecisal);
            mmagvti = mag(vintecisal);
            fprintf(1,'rintecisal hil %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n', ...
                rintecisal(1)*0.001,rintecisal(2)*0.001,rintecisal(3)*0.001,mmagrti*0.001,vintecisal(1)*0.001,vintecisal(2)*0.001,vintecisal(3)*0.001,mmagvti*0.001 );

            fprintf(1,'rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n', ...
                rhill2(1),rhill2(2),rhill2(3),vhill2(1),vhill2(2),vhill2(3) );
        end

    end % for ktr looping through each of 6 components

    fclose (outfilehill);

   

