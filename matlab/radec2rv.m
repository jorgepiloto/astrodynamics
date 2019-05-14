% ------------------------------------------------------------------------------
%
%                           function radec2rv
%
%  this function converts the right ascension and declination values with
%    position and velocity vectors of a satellite. uses velocity vector to
%    find the solution of singular cases.
%
%  author        : david vallado                  719-573-2600   25 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    rr          - radius of the satellite        er
%    rtasc       - right ascension                rad
%    decl        - declination                    rad
%    drr         - radius of the satellite rate   er/tu
%    drtasc      - right ascension rate           rad/tu
%    ddecl       - declination rate               rad/tu
%
%  outputs       :
%    r           -  position vector            er
%    v           -  velocity vector            er/tu
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
% [r,v] = radec2rv( rr,rtasc,decl,drr,drtasc,ddecl );
% ------------------------------------------------------------------------------

function [r,v] = radec2rv( rr,rtasc,decl,drr,drtasc,ddecl );

        % -------------------------  implementation   -----------------
        small        = 0.00000001;

        r(1)= rr*cos(decl)*cos(rtasc);
        r(2)= rr*cos(decl)*sin(rtasc);
        r(3)= rr*sin(decl);
        r = r';
        
        v(1)= drr*cos(decl)*cos(rtasc) - rr*sin(decl)*cos(rtasc)*ddecl ...
                 - rr*cos(decl)*sin(rtasc)*drtasc;
        v(2)= drr*cos(decl)*sin(rtasc) - rr*sin(decl)*sin(rtasc)*ddecl ...
                 + rr*cos(decl)*cos(rtasc)*drtasc;
        v(3)= drr*sin(decl) + rr*cos(decl)*ddecl;
        v = v';
