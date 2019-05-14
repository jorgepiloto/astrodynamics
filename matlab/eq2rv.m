% ------------------------------------------------------------------------------
%
%                           function eq2rv
%
%  this function finds the classical orbital elements given the equinoctial
%    elements.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%    vallado     - fix elliptical equatorial orbits case         19 oct 2002
%    vallado     - add constant file use                         29 jun 2003
%
%  inputs          description                    range / units
%    a           - semimajor axis                 km
%    af          - component of ecc vector
%    ag          - component of ecc vector
%    chi         - component of node vector in eqw
%    psi         - component of node vector in eqw
%    meanlon     - mean longitude                 rad
%
%  outputs       :
%    r           - position vector                km
%    v           - velocity vector                km/s
%
%  locals        :
%    n           - mean motion                    rad
%    temp        - temporary variable
%    p           - semilatus rectum               km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omega       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    m           - mean anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
%
%  coupling      :
%
%  references    :
%    vallado 2013:108
%
% [r, v] = eq2rv( a, af, ag, chi, psi, meanlon, fr);
% ------------------------------------------------------------------------------

function [r, v] = eq2rv( a, af, ag, chi, psi, meanlon, fr)

        % -------------------------  implementation   -----------------
        constmath;
        constastro;

        arglat  = 999999.1;
        lonper  = 999999.1;
        truelon = 999999.1;

        % ---- if n is input ----
        %a = (mu/n^2)^(1.0/3.0);

        ecc = sqrt (af^2 + ag^2);
        p = a * (1.0 - ecc*ecc);
        incl = pi*((1.0 - fr)*0.5) + 2.0*fr*atan( sqrt(chi^2 + psi^2) );
        omega = atan2( chi, psi);
        argp = atan2( ag,af ) - fr*atan2( chi,psi );

        if ( ecc < small )
            % ----------------  circular equatorial  ------------------
            if (incl<small) || ( abs(incl-pi)< small )
                argp = 0.0;
                omega= 0.0;
%                truelon = nu;
              else
                % --------------  circular inclined  ------------------
                argp= 0.0;
%                arglat = nu;
            end
          else
            % ---------------  elliptical equatorial  -----------------
            if ( ( incl<small) || (abs(incl-pi)<small) )
%                argp = lonper;
                omega= 0.0;
            end
        end
        
        m = meanlon - fr*omega - argp;
        m = rem (m + twopi,twopi);

        [e0, nu] = newtonm ( ecc, m );

        % ----------  fix for elliptical equatorial orbits ------------
        if ( ecc < small )
           % ----------------  circular equatorial  ------------------
           if (incl<small) || ( abs(incl-pi)< small )
               argp    = undefined;
               omega   = undefined;
               truelon = nu;
             else
               % --------------  circular inclined  ------------------
               argp  = undefined;
               arglat= nu;
           end
           nu   = undefined;
         else
           % ---------------  elliptical equatorial  -----------------
           if ( ( incl < small) || (abs(incl-pi) < small) )
               lonper = argp;
               argp    = undefined;
               omega   = undefined;
           end
         end

        % -------- now convert back to position and velocity vectors
        [r, v] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);

