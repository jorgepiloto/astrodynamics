%     -----------------------------------------------------------------
%
%                              Ex3_5.m
%
%  this file demonstrates example 3-5.
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

       % -------- gstime       - greenwich sidereal time
       year = 1992;
       mon = 8;
       day = 20;
       hr = 12;
       min = 14;
       sec = 0.00;
       [jdut1, jdut1frac] = jday( year,mon,day,hr,min,sec );

       fprintf(1,'jd %18.10f  %18.10f   %18.10f \n',jdut1, jdut1frac, jdut1+jdut1frac);

       % ---- gst = gstime(jdut1);
       tut1= ( jdut1+jdut1frac - 2451545.0 ) / 36525.0

       temp = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1  ...
              + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841

       % 360/86400 = 1/240, to deg, to rad
%       temp = rem( temp*rad/240.0,twopi )
       temp = rem( temp,86400 )

       temp = temp / 240.0

       % ------------------------ check quadrants --------------------
       if ( temp < 0.0 )
           temp = temp + 360;
         end

       gst = temp


       fprintf(1,'gst = %16.12f\n', gst );

       % -------- lstime       - local sidereal time
       lon   = -104.000/rad;
       fprintf(1,'input lon = %11.7f ',lon*rad );

       [lst,gst] = lstime ( lon, jdut1+jdut1frac );

       fprintf(1,'lst %14.10f gst %14.10f deg\n',lst*rad,gst*rad );

