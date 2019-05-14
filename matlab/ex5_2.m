%     -----------------------------------------------------------------
%
%                              Ex5_2.m
%
%  this file demonstrates example 5-2 and 5-4 for sun and moon rise/set.
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

        % --------  sun         - sun rise set
        [jd,jdfrac] = jday( 1996, 3, 23, 0, 0, 0.00 );
        latgd = 40.0/rad;
        lon = 0.00 / rad;
        whichkind = 's';

        [utsunrise,utsunset,error] = sunriset(jd+jdfrac,latgd,lon,whichkind);

        fprintf(1,'sun sunrise %14.4f  %14.4f  sunset %14.4f %14.4f \n',utsunrise, (utsunrise-floor(utsunrise))*60, utsunset, (utsunset-floor(utsunset))*60);


        [jd,jdfrac] = jday( 2011, 6, 25, 0, 0, 0.00 );
        latgd = 40.9/rad;
        lon = -74.3 / rad;
        whichkind = 's';

        [utsunrise,utsunset,error] = sunriset(jd+jdfrac,latgd,lon,whichkind);

        fprintf(1,'sun sunrise %14.4f  %14.4f  sunset %14.4f %14.4f \n',utsunrise, (utsunrise-floor(utsunrise))*60, utsunset, (utsunset-floor(utsunset))*60);

