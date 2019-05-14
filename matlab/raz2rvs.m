% ------------------------------------------------------------------------------
%
%                           function raz2rvs
%
%  this function converts range, azimuth, and elevation values with slant
%    range and velocity vectors for a satellite from a radar site in the
%    topocentric horizon (sez) system.
%
%  author        : david vallado                  719-573-2600   10 jun 2002
%
%  revisions
%    vallado     - del unnecessary code                          26 aug 2002
%
%  inputs          description                    range / units
%    rho         - satellite range from site      km
%    az          - azimuth                        0.0 to 2pi rad
%    el          - elevation                      -pi/2 to pi/2 rad
%    drho        - range rate                     km / s
%    daz         - azimuth rate                   rad / s
%    del         - elevation rate                 rad / s
%
%  outputs       :
%    rhovec      - sez satellite range vector     km
%    drhovec     - sez satellite velocity vector  km / s
%
%  locals        :
%    sinel       - variable for sin( el )
%    cosel       - variable for cos( el )
%    sinaz       - variable for sin( az )
%    cosaz       - variable for cos( az )
%    temp        -
%    temp1       -
%
%  coupling      :
%    mag         - magnitude of a vector
%
%  references    :
%    vallado       2001, 250-251, eq 4-4, eq 4-5
%
% [rhosez,drhosez] =  raz2rvs ( rho,az,el,drho,daz,del );
% ------------------------------------------------------------------------------

function [rhosez,drhosez] =  raz2rvs ( rho,az,el,drho,daz,del );

        % ----------------------- initialize values -------------------
        sinel= sin(el);
        cosel= cos(el);
        sinaz= sin(az);
        cosaz= cos(az);

        % ------------------- form sez range vector -------------------
        rhosez(1) = -rho*cosel*cosaz;
        rhosez(2) =  rho*cosel*sinaz;
        rhosez(3) =  rho*sinel;

        % ----------------- form sez velocity vector ------------------
        drhosez(1) = -drho*cosel*cosaz + rhosez(3)*del*cosaz + rhosez(2)*daz;
        drhosez(2) =  drho*cosel*sinaz - rhosez(3)*del*sinaz - rhosez(1)*daz;
        drhosez(3) =  drho*sinel       + rho*del*cosel;

