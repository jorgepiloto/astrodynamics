% ----------------------------------------------------------------------------
%
%                           function rv2eq.m
%
%  this function transforms a position and velocity vector into the flight
%    elements - latgc, lon, fpa, az, position and velocity magnitude.
%
%  author        : david vallado                  719-573-2600    7 jun 2002
%
%  revisions
%    vallado     - fix special orbit types (ee)                   5 sep 2002
%    vallado     - add constant file use                         29 jun 2003
%
%  inputs          description                    range / units
%    r           - eci position vector            km
%    v           - eci velocity vector            km/s
%
%  outputs       :
%    n           - mean motion                    rad
%    a           - semi major axis                km
%    af          - component of ecc vector
%    ag          - component of ecc vector
%    chi         - component of node vector in eqw
%    psi         - component of node vector in eqw
%    meanlon     - mean longitude                 rad
%    truelon     - true longitude                 rad
%
%  locals        :
%    none        -
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2013, 108
%    chobotov            30
%
% [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq ( r,v )
% ----------------------------------------------------------------------------

function [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq ( r,v )

        constmath;
        constastro;

        % -------- convert to classical elements ----------------------
        [p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r, v);
         
        % -------- setup retrograde factor ----------------------------
        fr = 1.0;
        % ---------- set this so it for orbits over 90 deg !! ---------
         if incl > pi*0.5
             fr = -1.0;
         end

        if ( ecc < small )
            % ----------------  circular equatorial  ------------------
            if (incl<small) || ( abs(incl-pi)< small )
                argp = 0.0;
                omega= 0.0;
%                nu   = truelon;
              else
                % --------------  circular inclined  ------------------
                argp= 0.0;
%                nu  = arglat;
            end
          else
            % ---------------  elliptical equatorial  -----------------
            if ( ( incl<small) || (abs(incl-pi)<small) )
                argp = lonper;
                omega= 0.0;
            end
        end

        af = ecc * cos(fr*omega + argp);
        ag = ecc * sin(fr*omega + argp);

        if (fr > 0 )
            chi = tan(incl*0.5) * sin(omega);
            psi = tan(incl*0.5) * cos(omega);
          else
            chi = cot(incl*0.5) * sin(omega);
            psi = cot(incl*0.5) * cos(omega);
        end;

        n  = sqrt(mu/(a*a*a));  

        meanlonM = fr*omega + argp + m; 
        meanlonM = rem(meanlonM,2.0*pi);

        meanlonNu = fr*omega + argp + nu;  
        meanlonNu = rem(meanlonNu,2.0*pi);
