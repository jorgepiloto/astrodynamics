%     -----------------------------------------------------------------
%
%                              Ex3_6.m
%
%  this file demonstrates example 3-6.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2007
%                            by david vallado
%
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

       % -------- gstime       - greenwich sidereal time
       year = 1992;
       mon = 8;
       day = 20;
       hr = 12;
       min = 14;
       sec = 0.00;
       [jdut1,jdut1frac] = jday( year,mon,day,hr,min,sec);

       fprintf(1,'jd %18.10f  %18.10f   %18.10f \n',jdut1, jdut1frac, jdut1+jdut1frac);

       fprintf(1,'\n--------gstime0 test \n' );
       fprintf(1,'input year %4i \n',year);

       gst0 = gstime0(year);

       fprintf(1,'gst0 = %14.8f deg\n', gst0*rad );


       % -------- lstime       - local sidereal time
       lon   = -104.000/rad;
       fprintf(1,'input lon = %11.7f ',lon*rad );

       [lst,gst] = lstime ( lon, jdut1+jdut1frac );

       fprintf(1,'gst = %14.8f\n', gst*rad );
       fprintf(1,'lst %11.7f gst %11.7f deg\n',lst*rad,gst*rad );





