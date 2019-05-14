% ------------------------------------------------------------------------------
%
%                           function rv2razel
%
%  this function converts geocentric equatorial (eci) position and velocity
%    vectors into range, azimuth, elevation, and rates.  notice the value
%    of small as it can affect the rate term calculations. the solution uses
%    the velocity vector to find the singular cases. also, the elevation and
%    azimuth rate terms are not observable unless the acceleration vector is
%    available.
%
%  author        : david vallado                  719-573-2600   22 jun 2002
%
%  revisions
%    vallado     - add terms for ast calculation                 30 sep 2002
%    vallado     - update for site fixes                          2 feb 2004 
%
%  inputs          description                    range / units
%    reci        - eci position vector            km
%    veci        - eci velocity vector            km/s
%    rs          - eci site position vector       km
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
%    rho         - satellite range from site      km
%    az          - azimuth                        0.0 to 2pi rad
%    el          - elevation                      -pi/2 to pi/2 rad
%    drho        - range rate                     km/s
%    daz         - azimuth rate                   rad / s
%    del         - elevation rate                 rad / s
%
%  locals        :
%    rhoveci     - eci range vector from site     km
%    drhoveci    - eci velocity vector from site  km / s
%    rhosez      - sez range vector from site     km
%    drhosez     - sez velocity vector from site  km
%    wcrossr     - cross product result           km / s
%    earthrate   - eci earth's rotation rate vec  rad / s
%    tempvec     - temporary vector
%    temp        - temporary real*8 value
%    temp1       - temporary real*8 value
%    i           - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    rot3        - rotation about the 3rd axis
%    rot2        - rotation about the 2nd axis
%
%  references    :
%    vallado       2007, 268-269, alg 27
%
% [rho,az,el,drho,daz,del] = rv2razel ( reci,veci, latgd,lon,alt,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );
% ------------------------------------------------------------------------------

function [rho,az,el,drho,daz,del] = rv2razel ( reci,veci, latgd,lon,alt,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );

    halfpi = pi*0.5;
    small  = 0.00000001;

    % --------------------- implementation ------------------------
    % ----------------- get site vector in ecef -------------------
    [rs,vs] = site ( latgd,lon,alt );

    % -------------------- convert eci to ecef --------------------
    a = [0;0;0];
    [recef,vecef,aecef] = eci2ecef(reci,veci,a,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);
    % simplified - just use sidereal time rotation
    % thetasa= 7.29211514670698e-05 * (1.0  - 0.0/86400.0 );
    % omegaearth = [0; 0; thetasa;];
    % [deltapsi,trueeps,meaneps,omega,nut] = nutation(ttt,ddpsi,ddeps);
    % [st,stdot] = sidereal(jdut1,deltapsi,meaneps,omega,0,0 );
    %  recef  = st'*reci;
    %  vecef  = st'*veci - cross( omegaearth,recef );


    % ------- find ecef range vector from site to satellite -------
    rhoecef  = recef - rs;
    drhoecef = vecef;
    rho      = mag(rhoecef);

    % ------------- convert to sez for calculations ---------------
    [tempvec]= rot3( rhoecef, lon          );
    [rhosez ]= rot2( tempvec, halfpi-latgd );

    [tempvec]= rot3( drhoecef, lon         );
    [drhosez]= rot2( tempvec,  halfpi-latgd);

    % ------------- calculate azimuth and elevation ---------------
    temp= sqrt( rhosez(1)*rhosez(1) + rhosez(2)*rhosez(2) );
    if ( ( temp < small ) )           % directly over the north pole
        el= sign(rhosez(3))*halfpi;   % +- 90 deg
    else
        magrhosez = mag(rhosez);
        el= asin( rhosez(3) / magrhosez );
    end

    if ( temp < small )
        az = atan2( drhosez(2), -drhosez(1) );
    else
        az= atan2( rhosez(2)/temp, -rhosez(1)/temp );
    end

    % ------ calculate range, azimuth and elevation rates ---------
    drho= dot(rhosez,drhosez)/rho;
    if ( abs( temp*temp ) > small )
        daz= ( drhosez(1)*rhosez(2) - drhosez(2)*rhosez(1) ) / ( temp*temp );
    else
        daz= 0.0;
    end

    if ( abs( temp ) > small )
        del= ( drhosez(3) - drho*sin( el ) ) / temp;
    else
        del= 0.0;
    end

