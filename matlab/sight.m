% ------------------------------------------------------------------------------
%
%                           function sight
%
%  this function takes the position vectors of two satellites and determines
%    if there is line-of-sight between the two satellites.  an oblate earth
%    with radius of 1 er is assumed.  the process forms the equation of
%    a line between the two vectors.  differentiating and setting to zero finds
%    the minimum value, and when plugged back into the original line equation,
%    gives the minimum distance.  the parameter tmin is allowed to range from
%    0.0  to 1.0 .  scale the k-component to account for oblate earth because it's
%    the only qunatity that changes.
%
%  author        : david vallado                  719-573-2600   31 oct 2003
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    r1          - position vector of the 1st sat km
%    r2          - position vector of the 2nd sat km
%    whichkind   - spherical or ellipsoidal earth 's', 'e'*default
%
%  outputs       :
%    los         - line of sight                  'yes','no '
%
%  locals        :
%    tr1         - scaled r1 vector
%    tr2         - scaled r2 vector
%    adotb       - dot product of a dot b
%    tmin        - minimum value of t from a to b
%    distsqrd    - min distance squared to earth
%    asqrd       - magnitude of a squared
%    bsqrd       - magnitude of b squared
%
%  coupling:
%
%  references    :
%    vallado       2001, 291-295, alg 35, ex 5-3
%
% [los] = sight ( r1,r2, whichkind );
% ------------------------------------------------------------------------------

function [los] = sight ( r1,r2, whichkind );

        eesqrd     =     0.006694385000;     % eccentricity of earth sqrd
        re = 6378.137; % km
        % -------------------------  implementation   -----------------
        for i=1 : 3
            tr1(i)= r1(i);
            tr2(i)= r2(i);
          end
        magr1 = mag(tr1);
        magr2 = mag(tr2);

        % --------------------- scale z component ---------------------
        if ( whichkind == 'e' )
            temp= 1.0 /sqrt(1.0 -eesqrd);
          else
            temp= 1.0;
          end
        tr1(3)= tr1(3)*temp;
        tr2(3)= tr2(3)*temp;
        bsqrd= magr2*magr2;
        asqrd= magr1*magr1;
        adotb= dot( tr1,tr2 );

        % ---------------------- find tmin ----------------------------
        distsqrd= 0.0;
        if ( abs(asqrd + bsqrd - 2.0 *adotb) < 0.0001  )
            tmin= 0.0;
          else
            tmin = ( asqrd - adotb ) / ( asqrd + bsqrd - 2.0 *adotb );
          end
tmin
        % ----------------------- check los ---------------------------
        if ( (tmin < 0.0 ) | (tmin > 1.0 ) )
            los= 'yes';
          else
            distsqrd= ( (1.0 -tmin)*asqrd + adotb*tmin )/re^2;
            if ( distsqrd > 1.0 )
                los= 'yes';
              else
                los= 'no ';
              end
          end
distsqrd

