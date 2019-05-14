% ------------------------------------------------------------------------------
%
%                           function rv2radec
%
%  this function converts the right ascension and declination values with
%    position and velocity vectors of a satellite. uses velocity vector to
%    find the solution of singular cases.
%
%  author        : david vallado                  719-573-2600   25 jun 2002
%
%  revisions
%    vallado     - fix rtasc tests                               29 sep 2002
%
%  inputs          description                    range / units
%    r           -  position vector               km
%    v           -  velocity vector               km/s
%
%  outputs       :
%    rr          - radius of the satellite        km
%    rtasc       - right ascension                rad
%    decl        - declination                    rad
%    drr         - radius of the satellite rate   km/s
%    drtasc      - right ascension rate           rad/s
%    ddecl       - declination rate               rad/s
%
%  locals        :
%    temp        - temporary position vector
%    temp1       - temporary variable
%
%  coupling      :
%    none
%
%  references    :
%    vallado       2001, 246-248, alg 25
%
% [rr,rtasc,decl,drr,drtasc,ddecl] = rv2radec( r,v );
% ------------------------------------------------------------------------------

function [rr,rtasc,decl,drr,drtasc,ddecl] = rv2radec( r,v );

        small  = 0.00000001;

        % -------------------------  implementation   -------------------------
        % ------------- calculate angles and rates ----------------
        rr= mag(r);
        temp= sqrt( r(1)*r(1) + r(2)*r(2) );
        if ( temp < small )
            rtasc= atan2( v(2) , v(1) );
          else
            rtasc= atan2( r(2) , r(1) );
          end
        decl= asin( r(3)/rr );

        temp1= -r(2)*r(2) - r(1)*r(1);  % different now
        drr= dot(r,v)/rr;
        if ( abs(temp1) > small )
            drtasc= ( v(1)*r(2) - v(2)*r(1) ) / temp1;
          else
            drtasc= 0.0;
          end
        if ( abs( temp ) > small )
            ddecl= ( v(3) - drr*sin( decl ) ) / temp;
          else
            ddecl= 0.0;
          end

