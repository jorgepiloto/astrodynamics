% ----------------------------------------------------------------------------
%
%                           function adbar2rv.m
%
%  this function transforms the adbarv elements (rtasc, decl, fpa, azimuth,
%    position and velocity magnitude) into eci position and velocity vectors.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    rmag        - eci position vector magnitude  km
%    vmag        - eci velocity vector magnitude  km/sec
%    rtasc       - right ascension of sateillite  rad
%    decl        - declination of satellite       rad
%    fpav        - sat flight path angle from vertrad
%    az          - sat flight path azimuth        rad
%
%  outputs       :
%    r           - eci position vector            km
%    v           - eci velocity vector            km/s
%
%  locals        :
%    none        -
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2001, xx
%    chobotov            70
%
% [r,v] = adbar2rv ( rmag,vmag,rtasc,decl,fpav,az );
% ----------------------------------------------------------------------------

function [r,v] = adbar2rv ( rmag,vmag,rtasc,decl,fpav,az );

        % -------- form position vector
        r(1)= rmag*cos(decl)*cos(rtasc);
        r(2)= rmag*cos(decl)*sin(rtasc);
        r(3)= rmag*sin(decl);

        % -------- form velocity vector
        v(1)= vmag*( cos(rtasc)*(-cos(az)*sin(fpav)*sin(decl) + ...
                     cos(fpav)*cos(decl)) - sin(az)*sin(fpav)*sin(rtasc) );
        v(2)= vmag*( sin(rtasc)*(-cos(az)*sin(fpav)*sin(decl) + ...
                     cos(fpav)*cos(decl)) + sin(az)*sin(fpav)*cos(rtasc) );
        v(3)= vmag*( cos(az)*cos(decl)*sin(fpav) + cos(fpav)*sin(decl) );

        r = r';
        v = v';

