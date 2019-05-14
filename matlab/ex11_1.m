%     -----------------------------------------------------------------
%
%                              Ex11_1.m
%
%  this file demonstrates example 11-1.
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
%             7 jun 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************

        constmath;

        % --------  satfov calculations
        incl = 40.0/rad;
        az   = 40.0/rad;
        slatgd = 50.0 /rad;
        slon   = 40.0 / rad;
        salt   = 800.0;  % km
        tfov   = 25.0 / rad;
        etactr = 0.0 / rad;
        [rhomin, rhomax] = ...
             satfov ( incl, az, slatgd, slon, salt,tfov,etactr );

        fprintf(1,'\noff axis test \n');
        % run twice with +- tfov angles to get the max and min from nadir
        % location
        incl = 40.0/rad;
        az   = 40.0/rad;
        slatgd = 50.0 /rad;
        slon   = 40.0 / rad;
        salt   = 800.0;  % km
        tfov   = 25.0 / rad;
        etactr = 40.0 / rad;
        [rhomin, rhomax] = ...
             satfov ( incl, az, slatgd, slon, salt,tfov,etactr );

        fprintf(1,'\noff axis test 2nd half \n');
        incl = 40.0/rad;
        az   = 40.0/rad;
        slatgd = 50.0 /rad;
        slon   = 40.0 / rad;
        salt   = 800.0;  % km
        tfov   = -25.0 / rad;
        etactr = 40.0 / rad;
        [rhomin, rhomax] = ...
             satfov ( incl, az, slatgd, slon, salt,tfov,etactr );
         
        fprintf(1,'\noff axis test circle calcs \n');
        incl = 40.0/rad;
        az   = 140.0/rad;
        slatgd = 50.0 /rad;
        slon   = 40.0 / rad;
        salt   = 800.0;  % km
        tfov   = 25.0 / rad;
        etactr = 27.5 / rad;
        [rhomin, rhomax] = ...
             satfov ( incl, az, slatgd, slon, salt,tfov,etactr );
         
        % ----- loop around the new circle with the sensor range ------
        for i= 0 : 18
            az= i*20.0 /rad;
            [tgtlat,tgtlon] = pathm( slatgd, slon, rhomin*0.5, az );
            if ( i == 0 )
                maxlat= tgtlat;
              end
            if ( i == 9 )
                minlat= tgtlat;
            end
           fprintf(1,'az %11.7f lat %11.7f lon %11.7f\n', az*rad, tgtlat*rad, tgtlon*rad );
            
        end

