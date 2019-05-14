% -----------------------------------------------------------------------------
%
%                           function sunriset
%
%  this function finds the universal time for sunrise and sunset given the
%    day and site location. use (- lon) for local time
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    jd          - julian date                    days from 4713 bc
%    latgd       - site latitude (south -)        -65ø to 65ø rad
%    lon         - site longitude (west -)        -2pi to 2pi rad
%    whichkind   - character for which rise/set   's' 'c' 'n' 'a'
%
%  outputs       :
%    utsunrise   - universal time of sunrise      hrs
%    utsunset    - universal time of sunset       hrs
%    error       - error parameter
%
%  locals        :
%    sunangle    - angle between the sun vector
%                  and a point on the earth     rad
%    jdtemp      - julian date for sunrise/set    days from 4713 bc
%    uttemp      - temporary ut time              days
%    tut1        - julian centuries from the
%                  jan 1, 2000 12 h epoch (ut1)
%    ra          - right ascension                rad
%    decl        - declination                    rad
%    meanlonsun  -                                rad
%    meananomalysun                               rad
%    lonecliptic - longitude of the ecliptic      rad
%    obliquity   - obliquity of the ecliptic      rad
%    gst         - for 0 h utc of each day        rad
%    lha         - local hour angle               rad
%    year        - year                           1900 .. 2100
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - hour                           0 .. 23
%    min         - minute                         0 .. 59
%    sec         - second                         0.0d0 .. 59.999d0
%    opt         - idx to for rise and set calc    1,2
%
%  coupling      :
%    invjulian   - finds the year day mon hr min sec from the julian date
%    julianday   - finds the julian date given year, mon day, hr, min, sec
%
%  references    :
%    vallado       2007, 283, Alg 30, Ex 5-2
%
% [utsunrise,utsunset,error] = sunriset(jd,latgd,lon,whichkind)
% -----------------------------------------------------------------------------

function [utsunrise,utsunset,error] = sunriset(jd,latgd,lon,whichkind)

        % ------------------------  implementation   ------------------
        twopi   =     2.0*pi;
        rad2deg =    180.0/pi;
        deg2rad =     pi/180.0;

        % -------------- make sure lon is within +- 180 deg -----------
        if ( lon > pi )
           lon = lon - twopi;
        end
        if ( lon < -pi )
           lon = lon + twopi;
        end

        if (whichkind == 's')
           sunangle = (90.0 + 50.0 / 60.0 ) * deg2rad;
        end
        if (whichkind == 'c')
           sunangle = 96.0 * deg2rad;
        end
        if (whichkind == 'n')
           sunangle = 102.0 * deg2rad;
        end
        if (whichkind == 'a')
           sunangle = 108.0 * deg2rad;
        end
        [year,month,day,hr,min,sec] = invjday(jd,0.0);

        for opt = 1:2,
           error   = 'ok';
            if ( opt == 1 )
               [jdtemp,jdtempf] = jday(year,month,day,6,0,0.0);
             else
               [jdtemp,jdtempf] = jday(year,month,day,18,0,0.0);
            end

            jdtemp = jdtemp +jdtempf - lon * rad2deg / 15.0 / 24.0;
jdtemp
lon * rad2deg / 15.0 / 24.0
            tut1 = (jdtemp - 2451545.0)/36525.0;
            meanlonsun = 280.4606184 + 36000.77005361 * tut1;
%            meanlonsun = 280.460 + 36000.770 * tut1;
            meananomalysun = 357.5277233 + 35999.05034 * tut1;
            meananomalysun = rem( meananomalysun * deg2rad,twopi );
            if ( meananomalysun < 0.0 )
               meananomalysun = meananomalysun + twopi;
            end
            lonecliptic = meanlonsun + 1.914666471 * sin(meananomalysun) ...
                        + 0.019994643 * sin(2.0 * meananomalysun);
            lonecliptic = rem( lonecliptic * deg2rad,twopi );
            if ( lonecliptic < 0.0 )
               lonecliptic = lonecliptic + twopi;
            end
            obliquity = 23.439291 - 0.0130042 * tut1;
            obliquity = obliquity * deg2rad;

        fprintf(1,'lonecl %11.7f tut1 %11.7f  obl %11.7f  \n', ...
                   lonecliptic*rad2deg,tut1,obliquity*rad2deg );

            ra   = atan( cos(obliquity) * tan(lonecliptic) );
            decl = asin( sin(obliquity) * sin(lonecliptic) );
            if ( ra < 0.0 )
               ra = ra + twopi;
            end
            if ( (lonecliptic > pi) & (ra < pi) )
               ra = ra + pi;
            end
            if ( (lonecliptic < pi) & (ra > pi) )
               ra = ra - pi;
            end
        fprintf(1,'mlonsun %11.7f meanansun %11.7f  eclon %11.7f  \n', ...
                   meanlonsun+1080.0,meananomalysun*rad2deg,lonecliptic*rad2deg );

        fprintf(1,'ra %11.7f decl %11.7f  \n', ra*rad2deg,decl*rad2deg );

            lha = (cos(sunangle) - sin(decl)*sin(latgd)) / (cos(decl)*cos(latgd) );
        %fprintf(1,'lha 1st %11.7f   \n', 90.0 + lha*rad2deg );
            if ( abs(lha) <= 1.0 )
               lha = acos( lha );
             else
               error = 'not ok';
       fprintf(1,'error \n');
            end
        fprintf(1,'lha 2nd  %11.7f   \n', lha*rad2deg );
            if ( error == 'ok' )
               if ( opt == 1 )
                  lha = twopi - lha;
               end
               
               gst = 1.75336855923327 + 628.331970688841 * tut1 ...
                   + 6.77071394490334e-06 * tut1 * tut1 ...
                   - 4.50876723431868e-10 * tut1 * tut1 * tut1;
               gst = rem( gst,twopi );
               if ( gst < 0.0 )
                  gst = gst + twopi;
               end
        fprintf(1,'lha %11.7f gst %11.7f  \n', lha*rad2deg,(gst-twopi)*rad2deg );
               uttemp = lha + ra  - gst ;  % rad       use   - lon for local time

        fprintf(1,'gst %11.7f uttemp %11.7f uttemp %11.7f   \n', gst*rad2deg,uttemp*rad2deg,(twopi+uttemp)*rad2deg );

               uttemp = uttemp * rad2deg / 15.0;  % hr
               uttemp = rem( uttemp,24.0 );

               if ( uttemp < 0.0 )

                  uttemp = uttemp + 24.0;
                  error  = 'day before';
               end
               if ( uttemp > 24.0 )
                  uttemp = uttemp - 24.0
                  error = 'day after';
               end
             else
               uttemp = 99.99;
            end
            if ( opt == 1 )
               utsunrise = uttemp;
             else
               utsunset = uttemp;
            end
        end

