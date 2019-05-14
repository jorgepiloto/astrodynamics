% ------------------------------------------------------------------------------
%
%                           procedure nodeonly
%
%  this procedure calculates the delta v's for a change in longitude of
%    ascending node only.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    deltaomega  - change in node                 rad
%    ecc         - ecc of first orbit
%    vinit       - initial velocity vector        er/tu
%    fpa         - flight path angle              rad
%    incl        - inclination                    rad
%
%
%  outputs       :
%    ifinal      - final inclination              rad
%    deltav      - change in velocity             er/tu
%
%  locals        :
%    vfinal      - final velocity vector          er/tu
%    arglat      - argument of latitude           rad
%    arglat1     - final argument of latitude     rad
%    nuinit      - initial true anomaly           rad
%    theta       -
%
%  coupling      :
%    asin      - arc sine function
%    acos      - arc cosine function
%
%  references    :
%    vallado       2007, 349, alg 40, ex 6-5
% function [ifinal,deltav ] = nodeonly(iinit,ecc,deltaomega,vinit,fpa,incl);
% ----------------------------------------------------------------------------- }

function [ifinal,deltav ] = nodeonly(iinit,ecc,deltaomega,vinit,fpa,incl);
    rad = 57.29577951308230;
    if ecc > 0.0000001
        % ------------------------- elliptical --------------------- }
        theta = atan( sin(iinit) * tan(deltaomega) );
        ifinal= asin( sin(theta) / sin(deltaomega) );
        deltav= 2.0 * vinit * cos(fpa) * sin(0.5 * theta);

        arglat = pi * 0.5; % set at 90 deg }
        arglat1= acos( cos(incl) * sin(incl) * (1.0-cos(deltaomega)) ...
            / sin(theta) );
    else
        % -------------------------- circular ---------------------- }
        ifinal = incl;
        theta = acos( cos(iinit) * cos(iinit) + sin(iinit) * sin(iinit) * cos(deltaomega) );
        deltav= 2.0 * vinit * sin(0.5 * theta);

        arglat = acos( tan(iinit) * (cos(deltaomega) - cos(theta)) ...
            / sin(theta) );
        arglat1= acos( cos(incl) * sin(incl) * (1.0 - cos(deltaomega)) ...
            / sin(theta) );
    end;

     fprintf(1,' theta   %11.7f deg \n',theta*rad );
     fprintf(1,' arglat   %11.7f  %11.7f  \n',arglat*rad, arglat1*rad );

 
