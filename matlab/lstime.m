% -----------------------------------------------------------------------------
%
%                           function lstime
%
%  this function finds the local sidereal time at a given location.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    lon         - site longitude (west -)        -2pi to 2pi rad
%    jd          - julian date                    days from 4713 bc
%
%  outputs       :
%    lst         - local sidereal time            0.0 to 2pi rad
%    gst         - greenwich sidereal time        0.0 to 2pi rad
%
%  locals        :
%    none.
%
%  coupling      :
%    gstime        finds the greenwich sidereal time
%
%  references    :
%    vallado       2007, 194, alg 15, ex 3-5
%
% [lst,gst] = lstime ( lon, jd );
% -----------------------------------------------------------------------------

function [lst,gst] = lstime ( lon, jd );

        twopi  = 2.0*pi;

        % ------------------------  implementation   ------------------
        [gst] = gstime( jd );
        lst = lon + gst;

        % ----------------------- check quadrants ---------------------
        lst = rem( lst,twopi );
        if ( lst < 0.0 )
            lst= lst + twopi;
          end

