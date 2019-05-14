%     -----------------------------------------------------------------
%
%                              Ex3_3.m
%
%  this file demonstrates example 3-3.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2007
%                            by david vallado
%
%     (h)               email dvallado@msn.com
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            13 feb 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************
        constmath;

        % --------  ijk2ll       - position to lat lon alt almanac (fastest)
        r1=[6524.834 6862.875 6448.296];
        [latgc,latgd,lon,hellp] = ijk2ll ( r1 );
        fprintf(1,'ijk2ll gc %14.7f gd %14.7f %14.7f%14.7f\n',latgc*rad,latgd*rad,lon*rad,hellp );

        % --------  ijk2llb      - position to lat lon alt borkowski
        [latgc,latgd,lon,hellp] = ijk2llb ( r1 );
        fprintf(1,'ijk2ll gc %14.7f gd %14.7f %14.7f%14.7f\n',latgc*rad,latgd*rad,lon*rad,hellp );


