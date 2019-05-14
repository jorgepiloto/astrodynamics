%
% ------------------------------------------------------------------------------
%
%                           function razel2rv
%
%  this function converts range, azimuth, and elevation and their rates to
%    the geocentric equatorial (eci) position and velocity vectors.
%
%  author        : david vallado                  719-573-2600   30 may 2002
%
%  revisions
%    vallado     - add terms for ast calculation                 30 sep 2002
%
%  inputs          description                    range / units
%    rho         - satellite range from site      km
%    az          - azimuth                        0.0 to 2pi rad
%    el          - elevation                      -pi/2 to pi/2 rad
%    drho        - range rate                     km/s
%    daz         - azimuth rate                   rad / s
%    del         - elevation rate                 rad / s
%    rs          - ecef site position vector      km
%    latgd       - geodetic latitude              -pi/2 to pi/2 rad
%    lon         - longitude of site              -2pi to 2pi rad
%    alt         - altitude                       km
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    terms       - number of terms for ast calculation 0,2
%
%  outputs       :
%    reci        - eci position vector            km
%    veci        - eci velocity vector            km/s
%
%  locals        :
%    rhoecef     - ecef range vector from site    km
%    drhoecef    - ecef velocity vector from site km/s
%    rhosez      - sez range vector from site     km
%    drhosez     - sez velocity vector from site  km
%    tempvec     - temporary vector
%
%  coupling      :
%    raz2rvs     - find r and v from site in topocentric horizon (sez) system
%
%  references    :
%    vallado       2001, 250-255, alg 27
%
% [reci,veci] = razel2rv ( rho,az,el,drho,daz,del,latgd,lon,alt,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );
% ------------------------------------------------------------------------------

function [reci,veci] = razel2rv ( rho,az,el,drho,daz,del,latgd,lon,alt,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );

        % -------------------------  implementation   -----------------
        constmath;

        % -----------  find sez range and velocity vectors ------------
        [rhosez,drhosez] = raz2rvs( rho,az,el,drho,daz,del );

        % -----------  perform sez to ijk (ecef) transformation -------
        [tempvec] = rot2( rhosez , latgd-halfpi );
        [rhoecef] = rot3( tempvec,-lon          );
        rhoecef = rhoecef';
        
        [tempvec] = rot2( drhosez, latgd-halfpi );
        [drhoecef]= rot3( tempvec,-lon          );
        drhoecef = drhoecef';

        % ----------  find ecef range and velocity vectors -------------
        [rs,vs] = site ( latgd,lon,alt );
        recef = rhoecef + rs;
        vecef = drhoecef;

        % -------- convert ecef to eci
        recef = recef;
        vecef = vecef;
        a     = [0;0;0];
        [reci,veci,aeci] = ecef2eci(recef,vecef,a,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );

