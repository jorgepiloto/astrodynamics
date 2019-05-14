%
% ----------------------------------------------------------------------------
%
%                           function rv2adbar.m
%
%  this function transforms a position and velocity vector into the adbarv
%    elements - rtasc, decl, fpav, azimuth, position and velocity magnitude.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    r           - eci position vector            km
%    v           - eci velocity vector            km/s
%
%  outputs       :
%    rmag        - eci position vector magnitude  km
%    vmag        - eci velocity vector magnitude  km/sec
%    rtasc       - right ascension of sateillite  rad
%    decl        - declination of satellite       rad
%    fpav        - sat flight path angle from vertrad
%    az          - sat flight path azimuth        rad
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
% [rmag,vmag,rtasc,decl,fpav,az] = rv2adbar ( r,v );
% ----------------------------------------------------------------------------

function [rmag,vmag,rtasc,decl,fpav,az] = rv2adbar ( r,v );

        twopi = 2.0*pi;
        small = 0.00000001;

        rmag = mag(r);

        vmag = mag(v);

        % ---------------- calculate rtasc and decl -------------------
        temp= sqrt( r(1)*r(1) + r(2)*r(2) );
        if ( temp < small )
            temp1= sqrt( v(1)*v(1) + v(2)*v(2) );
            if ( abs(temp1) > small )
                rtasc= atan2( v(2) , v(1) );
              else
                rtasc= 0.0;
              end
          else
            rtasc= atan2( r(2), r(1) );
          end
        decl= asin( r(3)/rmag );

        h    = cross(r,v);
        hmag = mag(h);
        rdotv= dot(r,v);
        fpav = atan2(hmag,rdotv);

        hcrossr = cross(h,r);
        az = atan2( r(1)*hcrossr(2) - r(2)*hcrossr(1), hcrossr(3)*rmag );

