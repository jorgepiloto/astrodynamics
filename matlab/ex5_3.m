%     -----------------------------------------------------------------
%
%                              Ex5_3.m
%
%  this file demonstrates example 5-3.
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

        % --------  moon         - analytical moon ephemeris
        jd = 2449470.5;
        [rmoon,rtasc,decl] = moon ( jd );
        fprintf(1,'moon rtasc %14.4f deg decl %14.4f deg\n',rtasc*rad,decl*rad );
        fprintf(1,'moon %14.7f%14.7f%14.7f er\n',rmoon );
        fprintf(1,'moon %14.4f%14.4f%14.4f km\n',rmoon*6378.137 );

