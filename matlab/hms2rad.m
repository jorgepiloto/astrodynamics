% -----------------------------------------------------------------------------
%
%                           function hms2rad
%
%  this function converts hours, minutes and seconds into radians.  notice
%    the conversion 0.2617 is simply the radian equivalent of 15 degrees.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    hr          - hours                          0 .. 24
%    min         - minutes                        0 .. 59
%    sec         - seconds                        0.0 .. 59.99
%
%  outputs       :
%    hms         - result                         rad
%
%  locals        :
%    temp        - conversion from hours to rad   0.261799
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2007, 204, alg 19 alg 20, ex 3-9
%
% [hms] = hms2rad( hr,min,sec );
% -----------------------------------------------------------------------------

function [hms] = hms2rad( hr,min,sec );

        % ------------------------  implementation   ------------------
        temp = 15.0 * pi/180.0;

        hms = ( hr + min/60.0 + sec/3600.0 )*temp;

