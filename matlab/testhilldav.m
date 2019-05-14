%
% test the hills conversion programs
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
    juststk  = 1;  % 0 is just run stk, no intermediate stuff
    ang_step = 0.00000001;
    
    % case where you read in the int and tgt ephemerides
    ropt = 'm'; % matlab or stk

    % clear the contents out between runs w deletes existing data in the file
    outfile = fopen('d:/STKFiles Educational Files/Hills/testhill.out', 'wt');  % write 
    fclose(outfile);
    outfile = fopen('d:/STKFiles Educational Files/Hills/testhill.out', 'at');  % append

    outfiledet = fopen(strcat('d:/STKFiles Educational Files/Hills/testhilldetc', int2str(casenumo),'n', int2str(casetest), '.out'), 'wt');      % rewrite each time
    outfiledethp = fopen(strcat('d:/STKFiles Educational Files/Hills/testhilldethpc', int2str(casenumo),'n', int2str(casetest), '.out'), 'wt');  % rewrite each time
    outfiledetH = fopen(strcat('d:/STKFiles Educational Files/Hills/testhilldetHc', int2str(casenumo),'n', int2str(casetest), '.out'), 'wt');      % rewrite each time
    outfiledethpH = fopen(strcat('d:/STKFiles Educational Files/Hills/testhilldethpHc', int2str(casenumo),'n', int2str(casetest), '.out'), 'wt');  % rewrite each time

    outfilehill = fopen(strcat('d:/STKFiles Educational Files/Hills/thillc', int2str(casenumo),'n', int2str(casetest), '.out'), 'wt');      % rewrite each time
    
    % ------------------------- target satellite info ------------------------------- }
    fprintf(1,'=================== be sure to have d:/STKFiles Educational Files/Hills/testhilla.sc open ================\n');

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

%---        [rintecix, vintecix]   = eqcm2eci( rtgteci,vtgteci,hrokm,hvokm );     % find int pos eqcm
%---        [rhillx, vhillx]      = eci2eqcm( rtgteci,vtgteci,rintecix,vintecix ); % find how much this would be for true hills
        [rintecix, vintecix]   = EQCM_to_ECI_RTN_sal( rtgteci, vtgteci, hrokm, hvokm );     % find int pos eqcm
        [rhillx, vhillx]      = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci, rintecix,vintecix, outfilehill ); % find how much this would be for true hills
        
        %        [rintecix,vintecix]   = eqcm2eci( rtgteci,vtgteci,hrokm,hvokm );     % find int pos eqcm
        %        [rhillx, vhillx]      = eci2eqcm( rtgteci,vtgteci,rintecix,vintecix ); % find how much this would be for true hills

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

    if (juststk ~= 0) 
    
     
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

%---        [rintecisal,vintecisal]   = eqcm2eci( rtgteci,vtgteci,hrokm,hvokm );              % find int pos eqcm
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

%---        [rintecisal,vintecisal]   = eqcm2eci( rtgteci,vtgteci,hrokm,hvokm );              % find int pos eqcm
%---        [rhill2, vhill2]          = eci2eqcm( rtgteci,vtgteci,rintecisal,vintecisal ); % find how much this would be for true hills
%---        [rintecisal1,vintecisal1] = eqcm2eci( rtgteci,vtgteci,rhill2,vhill2 );            % find int that would have found this new hills
%---        [rhillsal, vhillsal]      = eci2eqcm( rtgteci,vtgteci,rintecisal1,vintecisal1 );  %
%---        [rhillsal1, vhillsal1]    = eci2eqcm( rtgteci,vtgteci,rintecisal,vintecisal );    %
        [rintecisal,vintecisal]       = EQCM_to_ECI_RTN_sal( rtgteci,vtgteci,hrokm,hvokm );     % find int pos eqcm
        [rhill2, vhill2]              = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci,rintecisal,vintecisal, outfilehill ); % find how much this would be for true hills
        [rintecisal1,vintecisal1]     = EQCM_to_ECI_RTN_sal( rtgteci,vtgteci,rhill2,vhill2 );     % find int pos eqcm
        [rhillsal, vhillsal]          = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci,rintecisal1,vintecisal1, outfilehill ); % find how much this would be for true hills
        [rhillsal1, vhillsal1]        = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci,rintecisal,vintecisal, outfilehill ); % find how much this would be for true hills
        
        %            hro = rhill2*1000;  % m, used for hills propagtion
        %            hvo = vhill2*1000;

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
pause;
    end  % running just stk
   

    % -----------------------------------------------------------------------------------------------
    % -------------------------------------- now do it in STK ---------------------------------------
    % -----------------------------------------------------------------------------------------------
    fprintf(1,'\n\n-------------------------------------------------- now in STK -------------------------------------------- \n');
    ktr8str  = '10.0  0.01';
    ktr9str  = '100.0  0.01';
    ktr10str = '1000.0  0.01';
    ktr11str = '10.0  0.1';
    ktr12str = '100.0  0.1';
    ktr13str = '1000.0  0.1';
    ktr14str = '10.0  1.0';
    ktr15str = '100.0  1.0';
    ktr16str = '1000.0  1.0';
    first = 0;
    
    for ktr = casetest:casetest  % 16:16 do 8:13 for xls Case graphs, do 10:10 or 13:13 for hills plots and comparisons
        if ktr == 7
            x = 100000.0;  % in m
            y = 2348.0;  % in m
            z = -345.0;  % in m
            xd = 1;
            yd = 2.340;
            zd = 0.3829;
            currktrstr = 'misc';
        end
        if ktr == 8  % case for book
            x = 10.0;  % in m
            y = 10.0;%100.0;  % in m
            z = 10.0;%100.0;  % in m
            xd = 0.0100;%0.10;
            yd = 0.0100;%0.10;
            zd = 0.0100;%0.10;

%             x = 1000.0;  % in m
%             y = 1000.0;%100.0;  % in m
%             z = 1000.0;%100.0;  % in m
%             xd = 0.1000;%0.10;
%             yd = 0.1000;%0.10;
%             zd = 0.1000;%0.10;

            currktrstr = '10.0  0.01';
        end
        if ktr == 9  % case for book
            x = 100.0;  % in m
            y = 100.0;%100.0;  % in m
            z = 100.0;%100.0;  % in m
            xd = 0.0100;%0.10;
            yd = 0.0100;%0.10;
            zd = 0.0100;%0.10;
            currktrstr = '100.0  0.01';
        end
        if ktr == 10  % case for book
            x = 1000.0;  % in m
            y = 1000.0;%100.0;  % in m
            z = 1000.0;%100.0;  % in m
            xd = 0.0100;%0.10;
            yd = 0.0100;%0.10;
            zd = 0.0100;%0.10;
            currktrstr = '1000.0  0.01';
        end
        if ktr == 11  % case for book
            x = 10.0;  % in m
            y = 10.0;%100.0;  % in m
            z = 10.0;%100.0;  % in m
            xd = 0.100;%0.10;
            yd = 0.100;%0.10;
            zd = 0.100;%0.10;
            scalef = x*y*z/1000000.0;
            currktrstr = '10.0  0.1';
        end
        if ktr == 12  % case for book
            x = 100.0;  % in m
            y = 100.0;%100.0;  % in m
            z = 100.0;%100.0;  % in m
            xd = 0.100;%0.10;
            yd = 0.100;%0.10;
            zd = 0.100;%0.10;
            currktrstr = '100.0  0.1';
        end
        if ktr == 13  % case for book
            x = 1000.0;  % in m
            y = 1000.0;%100.0;  % in m
            z = 1000.0;%100.0;  % in m
            xd = 0.100;%0.10;
            yd = 0.100;%0.10;
            zd = 0.100;%0.10;
            currktrstr = '1000.0  0.1';
        end
        if ktr == 14  % case for book
            x = 10.0;  % in m
            y = 10.0;%100.0;  % in m
            z = 10.0;%100.0;  % in m
            xd = 1.00;%0.10;
            yd = 1.00;%0.10;
            zd = 1.00;%0.10;
            currktrstr = '10.0  1.0';
        end
        if ktr == 15  % case for book
            x =   100.0;  % in m
            y =   0.0;%100.0;  % in m
            z =   0.0;%100.0;  % in m
            xd = 0.00;%0.10;
            yd = 0.00;%0.10;
            zd = 0.00;%0.10;
            currktrstr = '100.0  1.0';
        end
        if ktr == 16  % case for book
            x = 10000.0;  % in m
            y = 10000.0;%100.0;  % in m
            z = 10000.0;%100.0;  % in m
            xd = 1.00;%0.10;
            yd = 1.00;%0.10;
            zd = 1.00;%0.10;
            currktrstr = '1000.0  1.0';
        end

        hro1(1)=  x/1000.0; % km
        hro1(2)=  y/1000.0;
        hro1(3)=  z/1000.0;
        hrokm = hro1'; % get in col format
        hvo1(1)= xd/1000.0; % km/s
        hvo1(2)= yd/1000.0;
        hvo1(3)= zd/1000.0;
        hvokm = hvo1'; % get in col format

%        [rintecisal,vintecisal]   = eqcm2eci( rtgteci,vtgteci,hrokm,hvokm );              % find int pos eqcm
%old way        [rintecisal,vintecisal]   = hill2eci( rtgteci,vtgteci,hrokm,hvokm );              % find int pos eqcm
 
        [rintecisal, vintecisal]   = EQCM_to_ECI_RTN_sal( rtgteci,vtgteci,hrokm,hvokm )     %  ,ang_step find int pos eqcm

%        pause;
%---        [rhill2, vhill2]          = eci2eqcm( rtgteci,vtgteci,rintecisal,vintecisal ); % find how much this would be for true hills
        [rhill2, vhill2]      = ECI_to_EQCM_RTN_sal( rtgteci,vtgteci,rintecisal,vintecisal, outfilehill  ); % ,ang_step find how much this would be for true hills

        % ------------ set the vectors for future use in program
        rinteci(1) = rintecisal(1);  % in km
        rinteci(2) = rintecisal(2);
        rinteci(3) = rintecisal(3);
        vinteci(1) = vintecisal(1);
        vinteci(2) = vintecisal(2);
        vinteci(3) = vintecisal(3);

        % set initial values to compare and make sure STK doesn't have an error
        rintecio = rinteci;  % in km 
        vintecio = vinteci;
        rtgtecio = rtgteci;
        vtgtecio = vtgteci;

%        if ktr >= 8  % cases for paper
            % get the vectors to STK, produce .e files, and then process the
            % differences for excel
            if first == 0
                fprintf(outfile,' ============================ Case for book %i \n',casenumo );
                mmagrti = mag(rinteci);
                mmagvti = mag(vinteci);
                fprintf(outfile,'rint      %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n', ...
                    rinteci(1),rinteci(2),rinteci(3),mmagrti,vinteci(1),vinteci(2),vinteci(3),mmagvti );
                mmagrti = mag(rtgteci);
                mmagvti = mag(vtgteci);
                fprintf(outfile,'rtgt      %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n', ...
                    rtgteci(1),rtgteci(2),rtgteci(3),mmagrti,vtgteci(1),vtgteci(2),vtgteci(3),mmagvti );
                fprintf(outfile,'rhill         %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f   \n', ...
                    rhill2(1),rhill2(2),rhill2(3),vhill2(1),vhill2(2),vhill2(3) );
                fprintf(outfiledet,' ============================ Case for book %i \n',casenumo );
                mmagrti = mag(rinteci);
                mmagvti = mag(vinteci);
                fprintf(outfiledet,'rint, ,     %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %15.8f %15.8f \n', ...
                    rinteci(1),rinteci(2),rinteci(3),vinteci(1),vinteci(2),vinteci(3),mmagrti,mmagvti );
                mmagrti = mag(rtgteci);
                mmagvti = mag(vtgteci);
                fprintf(outfiledet,'rtgt,  ,    %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %15.8f %15.8f \n', ...
                    rtgteci(1),rtgteci(2),rtgteci(3),vtgteci(1),vtgteci(2),vtgteci(3),mmagrti,mmagvti );
                fprintf(outfiledet,'rhill ,  ,      %20.13f, %20.13f, %20.13f,   %20.13f, %20.13f, %20.13f   \n', ...
                    rhill2(1),rhill2(2),rhill2(3),vhill2(1),vhill2(2),vhill2(3) );
                fprintf(outfiledethp,' ============================ Case for book %i \n',casenumo );
                mmagrti = mag(rinteci);
                mmagvti = mag(vinteci);
                fprintf(outfiledethp,'rint,  ,    %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %15.8f %15.8f \n', ...
                    rinteci(1),rinteci(2),rinteci(3),vinteci(1),vinteci(2),vinteci(3),mmagrti,mmagvti );
                mmagrti = mag(rtgteci);
                mmagvti = mag(vtgteci);
                fprintf(outfiledethp,'rtgt,  ,    %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %15.8f %15.8f \n', ...
                    rtgteci(1),rtgteci(2),rtgteci(3),vtgteci(1),vtgteci(2),vtgteci(3),mmagrti,mmagvti );
                fprintf(outfiledethp,'rhill ,  ,      %20.13f, %20.13f, %20.13f,   %20.13f, %20.13f, %20.13f   \n', ...
                    rhill2(1),rhill2(2),rhill2(3),vhill2(1),vhill2(2),vhill2(3) );
                first = 1;
            else
                fprintf(outfile,' \n');
                fprintf(outfiledet,' \n');
                fprintf(outfiledethp,' \n');             
            end
            
            % ---- Attach to running STK (stkeducationalfiles\testhilla.sc)
            app = actxGetRunningServer( 'STK.application' );
            root = app.get( 'Personality2' );  % general control (new interface)
            scen = root.CurrentScenario;

            % ---------------------------------- set these up one time ----------------------------------
            dtsec    =  180.0;   % step size in seconds
            numsteps = 1440; % number of points to consider, so 1440 at 3 min is 3 days
            % note that the .e files go for 4 days, or 5760 min

            %          fprintf(1,'start setting stk sats  - please wait... \n' );
            StartT = '1 jun 2011 00:00:00.000';
            StopT = '5 jun 2011 00:00:00.000';
            sat = scen.Children.Item( 'InterceptorHPOP' );  % Get current satellite
            %sat.SetPropagatorType( 'ePropagatorHPOP' );
            sat.Propagator.StartTime = StartT;
            sat.Propagator.StopTime = StopT;
            sat.Propagator.Step = dtsec;
            sat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', rinteci(1), rinteci(2), rinteci(3), vinteci(1), vinteci(2), vinteci(3) );
            sat.Propagator.Propagate;
            % Export ephemeris
            ephtool = sat.ExportTools.GetEphemerisStkExportTool;
            ephtool.CoordinateSystem = 'eStkEphemCoordinateSystemICRF';
            ephtool.Export(  'd:\STKFiles Educational Files\Hills\inthpop.txt' );

            sat = scen.Children.Item( 'Interceptor' );  % Get current satellite
            sat.SetPropagatorType( 'ePropagatorTwoBody' );  % already set in scenario
            sat.Propagator.StartTime = StartT;
            sat.Propagator.StopTime = StopT;
            sat.Propagator.Step = dtsec;
            sat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', rinteci(1), rinteci(2), rinteci(3), vinteci(1), vinteci(2), vinteci(3) );
            sat.Propagator.Propagate;
            % Export ephemeris
            ephtool = sat.ExportTools.GetEphemerisStkExportTool;
            ephtool.CoordinateSystem = 'eStkEphemCoordinateSystemICRF';
            ephtool.Export(  'd:\STKFiles Educational Files\Hills\int.txt' );

            sat = scen.Children.Item( 'TargetHPOP' );  % Get current satellite
            %sat.SetPropagatorType( 'ePropagatorHPOP' );
            sat.Propagator.StartTime = StartT;
            sat.Propagator.StopTime = StopT;
            sat.Propagator.Step = dtsec;
            sat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', rtgteci(1), rtgteci(2), rtgteci(3), vtgteci(1), vtgteci(2), vtgteci(3) );
            sat.Propagator.Propagate;
            % Export ephemeris
            ephtool = sat.ExportTools.GetEphemerisStkExportTool;
            ephtool.CoordinateSystem = 'eStkEphemCoordinateSystemICRF';
            ephtool.Export(  'd:\STKFiles Educational Files\Hills\tgthpop.txt' );

            sat = scen.Children.Item( 'Target' );  % Get current satellite
            %sat.SetPropagatorType( 'ePropagatorHPOP' );
            sat.Propagator.StartTime = StartT;
            sat.Propagator.StopTime = StopT;
            sat.Propagator.Step = dtsec;
            sat.Propagator.InitialState.Representation.AssignCartesian( 'eCoordinateSystemICRF', rtgteci(1), rtgteci(2), rtgteci(3), vtgteci(1), vtgteci(2), vtgteci(3) );
            sat.Propagator.Propagate;
            % Export ephemeris
            ephtool = sat.ExportTools.GetEphemerisStkExportTool;
            ephtool.CoordinateSystem = 'eStkEphemCoordinateSystemICRF';
            ephtool.Export(  'd:\STKFiles Educational Files\Hills\tgt.txt' );

            fprintf(1,'case %i, set maually \n', casenumo );
            %     pause;
            %          fprintf(1,'done with setting stk sats  \n' );

            % goto compiled program!!
    fprintf(1,'\n\n------------------------------- done in stk, go to compiled program? ---------------------------------------- \n');

    
    
    pause;        
            [res1] = readdote ( 'd:\STKFiles Educational Files\Hills\int.txt' );  % all these files are in m
            [res2] = readdote ( 'd:\STKFiles Educational Files\Hills\inthpop.txt' );
            [res3] = readdote ( 'd:\STKFiles Educational Files\Hills\tgt.txt' );
            [res4] = readdote ( 'd:\STKFiles Educational Files\Hills\tgthpop.txt' );

            % Check if stk mixes up state for low ecc, etc
%             if (mag(rintecio-res1.pos(:,1)'*0.001) > 0.001)
%                 fprintf(1,'error in rinteci, set maually %20.13f %20.13f %20.13f  \n',rintecio,res1.pos(:,1)*0.001 );
%                 pause;
%                 %    sat = scen.Children.Item( 'Interceptor' );  % Get current satellite
%                 %    sat.Propagator.Propagate;
%                 [res1] = readdote ( 'd:\STK Educational Files\Hills\int.txt' );  % all these files are in meters
%             end
%             if (mag(rintecio-res2.pos(:,1)'*0.001) > 0.001)
%                 fprintf(1,'error in H rinteci, set maually %20.13f %20.13f %20.13f \n',rintecio,res2.pos(:,1)*0.001 );
%             end
%             if (mag(rtgtecio-res3.pos(:,1)*0.001) > 0.001)
%                 fprintf(1,'error in rtgteci, set maually %20.13f %20.13f %20.13f \n',rtgtecio,res3.pos(:,1)*0.001 );
%             end
%             if (mag(rtgtecio-res4.pos(:,1)*0.001) > 0.001)
%                 fprintf(1,'error in H rtgteci, set maually %20.13f %20.13f %20.13f  \n',rtgtecio,res4.pos(:,1)*0.001 );
%             end
%             if (mag(vintecio-res1.vel(:,1)'*0.001) > 0.1)
%                 fprintf(1,'error in vinteci, set maually %20.13f %20.13f %20.13f \n',vintecio,res1.vel(:,1)*0.001 );
%                 pause;
%                 %    sat = scen.Children.Item( 'Interceptor' );  % Get current satellite
%                 %    sat.Propagator.Propagate;
%                 [res1] = readdote ( 'd:\STK Educational Files\Hills\int.txt' );  % all these files are in meters
%             end
%             if (mag(vintecio-res2.vel(:,1)'*0.001) > 0.1)
%                 fprintf(1,'error in H vinteci, set maually %20.13f %20.13f %20.13f \n',vintecio,res2.vel(:,1)*0.001 );
%             end
%             if (mag(vtgtecio-res3.vel(:,1)*0.001) > 0.1)
%                 fprintf(1,'error in vtgteci, set maually %20.13f %20.13f %20.13f \n',vtgtecio,res3.vel(:,1)*0.001 );
%             end
%             if (mag(vtgtecio-res4.vel(:,1)*0.001) > 0.1)
%                 fprintf(1,'error in H vtgteci, set maually %20.13f %20.13f %20.13f \n', vtgtecio,res4.vel(:,1)*0.001 );
%             end

            fprintf(1,'start stepping through times and finding diffs  \n' );

            % setup initial hills displacement from various sources (above)
            kk = 1;
 
            % use the linear approximation only
            hrost = hrokm;  % km
            hvost = hvokm;
            ri = res3.pos(:,kk);  % in m, and must be the target 3!!!!
% xxx m doesn't work
%            hrost = hro;  % m
%            hvost = hvo;

            % use the two-body values 
%            [hrost, hvost] = eci2eqcm( res3.pos(:,kk),res3.vel(:,kk),res1.pos(:,kk),res1.vel(:,kk) ); % two-body
%             hrost = hrost / 1000;
%             hvost = hvost / 1000;
% units must be kmxxx hpop
%            ri = res3.pos(:,kk);  % must be the target!!!!

            % use the hpop values 
%            [hrost, hvost] = eci2eqcm(res4.pos(:,kk),res4.vel(:,kk),res2.pos(:,kk),res2.vel(:,kk) );  
%             hrost = hrost / 1000;
%             hvost = hvost / 1000;
% units must be kmxxx hpop
%            ri = res4.pos(:,kk);  % must be the target!!!!

            magralt = (mag(ri)-6378137.0)*0.001   % calc and change to km
            
            % select various initial hill's estimates
            for kk = 1 : numsteps %  1440 % 1 is acually the first point at 0 time...
                %load data into x y z arrays
                rinteci2 = res1.pos(:,kk)/1000.0;  % two-body in m, to km
                vinteci2 = res1.vel(:,kk)/1000.0;
                ttime = res1.t(:,kk);  % sec

                rintecih = res2.pos(:,kk)/1000.0;  % hpop in m
                vintecih = res2.vel(:,kk)/1000.0;

                rtgteci2 = res3.pos(:,kk)/1000.0;  % two-body in m
                vtgteci2 = res3.vel(:,kk)/1000.0;

                rtgtecih = res4.pos(:,kk)/1000.0;  % hpop in m
                vtgtecih = res4.vel(:,kk)/1000.0;

                % do conversions at each step from the numerical ephemerides
%---                [rhill2body, vhill2body] = eci2eqcm( rtgteci2*0.001, vtgteci2*0.001, rinteci2*0.001, vinteci2*0.001 ); % in km 
%---                [rhillhpop, vhillhpop]   = eci2eqcm( rtgtecih*0.001, vtgtecih*0.001, rintecih*0.001, vintecih*0.001 );
                [rhill2body, vhill2body]      = ECI_to_EQCM_RTN_sal( rtgteci2, vtgteci2, rinteci2, vinteci2, outfilehill ); % ,ang_step find how much this would be for true hills
                [rhillhpop, vhillhpop]      = ECI_to_EQCM_RTN_sal( rtgtecih, vtgtecih, rintecih, vintecih, outfilehill ); % ,ang_step find how much this would be for true hills
%rhill2body
                
% xxx
                % if sal, then change 1 more asignments below and hrost above... 
%               [rhill2body, vhill2body] = ECI_to_EQCM_NTW_sal( rtgteci2,vtgteci2,rinteci2,vinteci2 ); % in m 
%               [rhillhpop, vhillhpop]   = ECI_to_EQCM_NTW_sal( rtgtecih,vtgtecih,rintecih,vintecih );
               rhill2body = rhill2body'; % in km, chg orientation
               vhill2body = vhill2body'; % in km
               rhillhpop = rhillhpop'; % in km
               vhillhpop = vhillhpop'; % in km

                % use Hills with various initial conditions
                [rhill, vhill] = hillsr( hrost,hvost, magralt, ttime);  % in km here, eqcm in m or [km], alt in km, time in sec (1 in steps)

                % store rhill for later plotting against numerical versions
                rhillarr(kk,1) = rhill(1);  % this is in km because of input to hillrs
                rhillarr(kk,2) = rhill(2);
                rhillarr(kk,3) = rhill(3);
                rhillhpoparr(kk,1) = rhillhpop(1); % in km
                rhillhpoparr(kk,2) = rhillhpop(2);
                rhillhpoparr(kk,3) = rhillhpop(3);
                rhill2barr(kk,1) = rhill2body(1);
                rhill2barr(kk,2) = rhill2body(2);
                rhill2barr(kk,3) = rhill2body(3);

                vhillarr(kk,1) = vhill(1);
                vhillarr(kk,2) = vhill(2);
                vhillarr(kk,3) = vhill(3);
                vhillhpoparr(kk,1) = vhillhpop(1);
                vhillhpoparr(kk,2) = vhillhpop(2);
                vhillhpoparr(kk,3) = vhillhpop(3);
                vhill2barr(kk,1) = vhill2body(1);
                vhill2barr(kk,2) = vhill2body(2);
                vhill2barr(kk,3) = vhill2body(3);

                dr  = rhill2body - rhill;    % error in hills approx, all in km
                dv  = vhill2body - vhill;
                drh = rhillhpop - rhill;
                dvh = vhillhpop  - vhill;
              
if (rem(kk, 50) == 0) || (kk == 1)                
%if (kk>798) && (kk <801)  
    fprintf( 1,'%6i,%12.3f  h2b %6.3f %6.3f %6.3f, h %6.3f %6.3f %6.3f, dr %6.3f %6.3f %6.3f, int %11.3f %11.3f %11.3f, tgt %11.3f %11.3f %11.3f \n', ...
          kk, ttime, rhill2body(1),rhill2body(2),rhill2body(3),rhill(1),rhill(2),rhill(3), dr(1),dr(2),dr(3), ...
          rinteci2(1),rinteci2(2),rinteci2(3), rtgteci2(1),rtgteci2(2),rtgteci2(3) );
%     dbstop in testhilldav at 771 
end
%     1,   0.000  h2b  0.100  0.000  0.000, h  0.100  0.000  0.000, int 6878237.000       0.000       0.000, tgt 6878137.000       0.000       0.000 
%     2, 180.000  h2b  0.102 -0.020  0.000, h  0.106 -0.001  0.000, int 6742198.558 1361223.730      23.758, tgt 6742094.602 1361223.334      23.758 
%     3, 360.000  h2b  0.108 -0.042 -0.000, h  0.124 -0.006  0.000, int 6339464.620 2668602.690      46.576, tgt 6339348.963 2668599.503      46.576 
%     4, 540.000  h2b  0.117 -0.067 -0.000, h  0.152 -0.021  0.000, int 5685966.422 3870422.267      67.552, tgt 5685831.863 3870411.426      67.551 
%     5, 720.000  h2b  0.130 -0.096 -0.000, h  0.190 -0.049  0.000, int 4807554.692 4919143.853      85.855, tgt 4807395.079 4919117.922      85.855 


% xxx
                dr = dr * 1000; % in m
                dv = dv * 1000; % in m
                drh = drh * 1000; % in m
                dvh = dvh * 1000; % in m

%                 drpp1 = rinteci2 - rintecih; % overall difference in vectors, 2body and hpop
%                 drpp2 = rtgteci2 - rtgtecih;
%                 dvpp1 = vinteci2 - vintecih; % overall difference in vectors, 2body and hpop
%                 dvpp2 = vtgteci2 - vtgtecih;
%                 drp(kk,3) = mag(drpp1);
%                 drp(kk,4) = mag(drpp2);
%                 dvp(kk,3) = mag(dvpp1);
%                 dvp(kk,4) = mag(dvpp2);

                dra(kk,1) = dr(1); % m
                dra(kk,2) = dr(2);
                dra(kk,3) = dr(3);
                dra(kk,4) = mag(dr);
                drha(kk,1) = drh(1);
                drha(kk,2) = drh(2);
                drha(kk,3) = drh(3);
                drha(kk,4) = mag(drh);
                dva(kk,1) = dv(1); % m/s
                dva(kk,2) = dv(2);
                dva(kk,3) = dv(3);
                dva(kk,4) = mag(dv);
                dvha(kk,1) = dvh(1);
                dvha(kk,2) = dvh(2);
                dvha(kk,3) = dvh(3);
                dvha(kk,4) = mag(dvh);
                
               fprintf(outfiledet,'diff 2body, %6i, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f \n', ...
                      ttime, dr(1),dr(2),dr(3),dv(1),dv(2),dv(3) );
               fprintf(outfiledethp,'diff hpop,  %6i, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f \n', ...
                      ttime, drh(1),drh(2),drh(3),dvh(1),dvh(2),dvh(3) );
               fprintf(outfiledetH,'Hill hill,  %6i, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f \n', ...
                      ttime, rhill(1),rhill(2),rhill(3),vhill(1),vhill(2),vhill(3) );
               fprintf(outfiledethpH,'Hill hpop,  %6i, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f, %20.13f \n', ...
                      ttime, rhillhpop(1),rhillhpop(2),rhillhpop(3),vhillhpop(1),vhillhpop(2),vhillhpop(3) );

             %   [rhill2body, vhill2body]      = ECI_to_EQCM_RTN_sal( rtgteci2, vtgteci2, rinteci2, vinteci2 ); % find how much this would be for true hills
             rx(1,1) = rhill2body(1,1); % keep in km
             rx(2,1) = rhill2body(1,2); 
             rx(3,1) = rhill2body(1,3); 
             vx(1,1) = rhill2body(1,1); 
             vx(2,1) = rhill2body(1,2); 
             vx(3,1) = rhill2body(1,3);              
             [rintecisl, vintecisl]   = EQCM_to_ECI_RTN_sal( rtgteci2, vtgteci2, rx, vx );     % ,ang_step find int pos eqcm
              fprintf(outfilehill,' %6i %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n', ...
                      ttime, rintecisl(1),rintecisl(2),rintecisl(3),vintecisl(1),vintecisl(2),vintecisl(3) );
                  
                  
            end   % looping through the hpop and 2body cases...

    end % for ktr looping through each of 6 tests

  
fclose( outfiledet );
fclose( outfiledethp );
fclose( outfiledetH );
fclose( outfiledethpH );

fclose( outfilehill );
    
%pause;


%         res1.t(1:numsteps) = res1.t(1:numsteps) / 86400.0;
%         %res1.t(numsteps)
% 
%         %    figure(1);
%         %    plot(res1.t(1:numsteps),dra(:,4),'b-','LineWidth',1.5); % magnitude
%         %    hold on;
% 
% 
%         if ktr == 8 || ktr == 1
%             color1 = 'b-';
%         end;
%         if ktr == 9 || ktr == 2
%             color1 = 'c-';
%         end;
%         if ktr == 10 || ktr == 3
%             color1 = 'r-';
%         end;
%         if ktr == 11
%             color1 = 'm-';
%         end;
%         if ktr == 12
%             color1 = 'g-';
%         end;
%         if ktr == 13
%             color1 = 'y-';
%         end;
%         %    plot(res1.t(1:numsteps),dra(:,1),'b-'); %
%         %    plot(res1.t(1:numsteps),dra(:,2),'b-.'); %
%         %    plot(res1.t(1:numsteps),dra(:,3),'b-.'); %
%         figure(1);
%         %plot(res1.t(1:numsteps),dra(:,4),color1,'LineWidth',1.5); % magnitude
%         %hold on;
%         nummpoints = 10;
%         mytemp = reshape( dra(:,4), numsteps/nummpoints, nummpoints );
%         mytime = reshape( res1.t(1:numsteps), numsteps/nummpoints, nummpoints );
%         errorbar( mean(mytime), mean(mytemp), std(mytemp), color1 )
%         legend( ktr8str, ktr9str, ktr10str, ktr11str, ktr12str, ktr13str, 'Location', 'NorthWest' )
% 
%         fprintf(1,'kep  case %6.2f %6.2f  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n', ...
%             x, xd, mean(dra(1:240,4)),mean(dra(241:480,4)),mean(dra(481:720,4)),mean(dra(721:960,4)),mean(dra(961:1200,4)),mean(dra(1201:1440,4)) );
%         fprintf(1,'kep  case %6.2f %6.2f  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n', ...
%             x, xd, std(dra(1:240,4)),std(dra(241:480,4)),std(dra(481:720,4)),std(dra(721:960,4)),std(dra(961:1200,4)),std(dra(1201:1440,4)) );
%         p = polyfit( mean(mytime), mean(mytemp),1); % find the polynomial fit
%         fprintf(outfile,'kep  %s  p = %f (t) + %f \n', currktrstr,p(1),p(2) )
%         if ktr == 8
%             % xlabel('time (days)');
%             % ylabel('error (m)');
%             xlabel(gcf,'String','time (days)','FontName','Times New Roman','FontWeight','normal','FontSize',9);
%             ylabel(gcf,'String','error (m)','FontName','Times New Roman','FontWeight','normal','FontSize',9);
%             %        set(gca,'YLim',[1.0, 1000.0]);
%             % set(gca,'YScale','log');
%             set(gca,'XMinorTick','on');
%             set(gca,'YMinorTick','on');
%             hold on;
%             drawnow;
%         end
%         %     % ---- for first graph
%         %     figure;
%         %     plot(res1.t(1:numsteps),dra(:,1),'r-'); %
%         %     hold on
%         %     plot(res1.t(1:numsteps),dra(:,2),'c-'); %
%         %     plot(res1.t(1:numsteps),dra(:,3),'m-'); %
%         %     plot(res1.t(1:numsteps),dra(:,4),'b-'); %
%         %     numpoints = 10;
%         %     mytemp = reshape( dra(:,4), numsteps/nummpoints, nummpoints );
%         %     mytime = reshape( res1.t(1:numsteps), numsteps/nummpoints, nummpoints );
%         %     errorbar( mean(mytime), mean(mytemp), std(mytemp), 'g-','LineWidth',1.5 )
%         %     legend( 'Radial', 'Along-track','Cross-track','Magnitude', 'Errorbars', 'Location', 'NorthWest' )
% 
%         % check the tgt 2 tgt differences
%         %mytemp = reshape( drp(:,4), numsteps/nummpoints, nummpoints );
%         %errorbar( mean(mytime), mean(mytemp), std(mytemp), color1 );
% 
% 
% 
%         %    plot(res1.t(1:numsteps),drha(:,1),'m-'); %
%         %    plot(res1.t(1:numsteps),drha(:,2),'m-.'); %
%         %    plot(res1.t(1:numsteps),drha(:,3),'m-.'); %
%         figure(2);
%         %plot(res1.t(1:numsteps),drha(:,4),color1,'LineWidth',1.5); % magnitude
%         %hold on;
%         mytemp = reshape( drha(:,4), numsteps/nummpoints, nummpoints );
%         mytime = reshape( res1.t(1:numsteps), numsteps/nummpoints, nummpoints );
%         errorbar( mean(mytime), mean(mytemp), std(mytemp), color1 )
%         legend( ktr8str, ktr9str, ktr10str, ktr11str, ktr12str, ktr13str, 'Location', 'NorthWest' )
% 
%         fprintf(1,'hpop case %6.2f %6.2f  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n', ...
%             x, xd, mean(drha(1:240,4)),mean(drha(241:480,4)),mean(drha(481:720,4)),mean(drha(721:960,4)),mean(drha(961:1200,4)),mean(drha(1201:1440,4)) );
%         fprintf(1,'hpop case %6.2f %6.2f  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n', ...
%             x, xd, std(drha(1:240,4)),std(drha(241:480,4)),std(drha(481:720,4)),std(drha(721:960,4)),std(drha(961:1200,4)),std(drha(1201:1440,4)) );
%         p = polyfit( mean(mytime), mean(mytemp),1);
%         fprintf(outfile,'hpop %s  p = %f (t) + %f \n', currktrstr,p(1),p(2) )
%         if ktr == 8
%             % xlabel('time (days)');
%             % ylabel('error (m)');
%             xlabel(gcf,'String','time (days)','FontName','Times New Roman','FontWeight','normal','FontSize',9);
%             ylabel(gcf,'String','error (m)','FontName','Times New Roman','FontWeight','normal','FontSize',9);
%             %        set(gca,'YLim',[1.0, 1000.0]);
%             % set(gca,'YScale','log');
%             set(gca,'XMinorTick','on');
%             set(gca,'YMinorTick','on');
%             hold on;
%             drawnow;
%         end
% 
%         if ktr == 888
%             figure;
%             plot(rhillhpoparr(:,2),rhillhpoparr(:,1),'r-');
%             hold on;
%             plot(rhillarr(:,2),rhillarr(:,1),'b-');
% 
%             figure;
%             plot(rhillhpoparr(:,2),rhillhpoparr(:,1),'r-');
%             hold on;
%             plot(rhillarr(:,2),rhillarr(:,1),'b-');
%             axis equal;
% 
%             figure;
%             plot(rhillhpoparr(:,2),rhillhpoparr(:,3),'r-');
%             hold on;
%             plot(rhillarr(:,2),rhillarr(:,3),'b-');
% 
%             figure;
%             plot(vhillhpoparr(:,2),vhillhpoparr(:,1),'r-');
%             hold on;
%             plot(vhillarr(:,2),vhillarr(:,1),'b-');
% 
%             figure;
%             plot(vhillhpoparr(:,2),vhillhpoparr(:,1),'r-');
%             hold on;
%             plot(vhillarr(:,2),vhillarr(:,1),'b-');
%             axis equal;
% 
%             figure;
%             plot(vhillhpoparr(:,2),vhillhpoparr(:,3),'r-');
%             hold on;
%             plot(vhillarr(:,2),vhillarr(:,3),'b-');
% 
%         end
% 


% old end of loop...
    
    fprintf(1,'finished with diffs \n');


    % draw error bars for a segment and at a point in time
    %errorbar(.25,mean(dra(1:240,4)),std(dra(1:240,4)),'b-','LineWidth',1.5); % magnitude

    %  mytemp = reshape( dra(:,4), 1440/20, 20 );
    %  mytime = reshape( res1.t(1:numsteps), 1440/20, 20 );
    %  errorbar( mean(mytime), mean(mytemp), std(mytemp), 'b-' )
    %  legend( 'plot1', 'plot2','plot3','plot34', 'Location', 'NorthWest' )


    fprintf(1,'results in d:_STK Educational FIles_Hills_testhill.out  \n' );
    fprintf(1,'results in d:_STK Educational FIles_Hills_testhilldet.out  \n' );
    fprintf(1,'should be done now. sort the outfile and plot in excel \n' );

    fclose (outfile);

 %   load handel
 %   sound(y,Fs)
    
    
    
    