% -----------------------------------------------------------------------------
%
%                           function gstime0
%
%  this function finds the greenwich sidereal time at the beginning of a year.
%    this formula is derived from the astronomical almanac and is good only
%    0 hr ut1, jan 1 of a year.
%
%  author        : david vallado                  719-573-2600   22 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    year        - year                           1998, 1999, etc.
%
%  outputs       :
%    gst  0      - greenwich sidereal time        0 to 2pi rad
%
%  locals        :
%    jd          - julian date                    days from 4713 bc
%    temp        - temporary variable for reals   rad
%    tut1        - julian centuries from the
%                  jan 1, 2000 12 h epoch (ut1)
%
%  coupling      :
%
%  references    :
%    vallado        2007, 195, Eq 3-46
%
% gst0 = gstime0(year);
% -----------------------------------------------------------------------------

function gst0 = gstime0(year);

        twopi = 2.0*pi;

        % ------------------------  implementation   ------------------
        jd = 367.0 * year  ...
             - floor( (7 * (year + floor( 10 / 12.0) ) ) * 0.25 )   ...
             + floor( 275 / 9.0 ) ...
             + 1721014.5;   
         
        tut1 = ( jd - 2451545.0 ) / 36525.0;

        temp = - 6.2e-6 * tut1 * tut1 * tut1   ...
             + 0.093104 * tut1 * tut1  ...
             + (876600.0  * 3600.0 + 8640184.812866 ) * tut1  ...
             + 67310.54841;

        % ------------------------ check quadrants --------------------
        temp = rem( temp,twopi );
        if ( temp < 0.0 )
           temp = temp + twopi;
        end

        gst0 = temp;

