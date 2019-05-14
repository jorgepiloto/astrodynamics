%     -----------------------------------------------------------------
%
%                              Ex3_9.m
%
%  this file demonstrates example 3-9.
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

       % -------- hms test
       deg = 15;
       min = 15;
       sec = 53.63;
       fprintf(1,'deg %4i ',deg);
       fprintf(1,'min %4i ',min);
       fprintf(1,'sec %8.6f ',sec);

        [hms] = hms2rad( deg,min,sec );

       fprintf(1,'hms %11.7f \n',hms);

       [hr,min,sec] = rad2hms( hms );

       fprintf(1,' hr min sec %4i  %4i  %8.6f \n',hr, min, sec);

