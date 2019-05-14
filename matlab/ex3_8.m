%     -----------------------------------------------------------------
%
%                              Ex3_8.m
%
%  this file demonstrates example 3-8.
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

       % -------- dms test
       deg = -35;
       min = -15;
       sec = -53.63;
       fprintf(1,'deg %4i ',deg);
       fprintf(1,'min %4i ',min);
       fprintf(1,'sec %8.6f \n',sec);

       [dms] = dms2rad( deg,min,sec );

       fprintf(1,'dms %11.7f \n',dms);


       [deg,min,sec] = rad2dms( dms );

       fprintf(1,' deg min sec %4i  %4i  %8.6f \n',deg, min, sec);


