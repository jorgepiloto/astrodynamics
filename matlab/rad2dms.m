% -----------------------------------------------------------------------------
%
%                           function rad2dms
%
%  this function converts radians to degrees, minutes and seconds.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    dms         - result                         rad
%
%  outputs       :
%    deg         - degrees                        0 .. 360
%    min         - minutes                        0 .. 59
%    sec         - seconds                        0.0 .. 59.99
%
%  locals        :
%    temp        - temporary variable
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2001, 199, alg 17 alg 18, ex 3-8
%
% [deg,min,sec] = rad2dms( dms );
% -----------------------------------------------------------------------------

function [deg,min,sec] = rad2dms( dms );

        rad2deg = 180.0/pi;

        % ------------------------  implementation   ------------------
        temp = dms * rad2deg;
        deg = fix( temp );
        min = fix( (temp - deg)*60.0 );
        sec = (temp - deg - min/60.0 ) * 3600.0;

