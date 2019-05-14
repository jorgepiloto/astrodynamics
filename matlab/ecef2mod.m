%
% ----------------------------------------------------------------------------
%
%                           function ecef2mod
%
%  this function transforms a vector from the earth fixed (itrf) frame, to
%    the mean of date (mod) frame.
%
%  author        : david vallado                  719-573-2600    4 jun 2002
%
%  revisions
%
%  inputs          description                    range / units
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    eqeterms    - terms for ast calculation      0,2
%    ddpsi       - delta psi correction to gcrf   rad
%    ddeps       - delta eps correction to gcrf   rad
%
%  outputs       :
%    rmod        - position vector mod            km
%    vmod        - velocity vector mod            km/s
%    amod        - acceleration vector mod        km/s2
%
%  locals        :
%    deltapsi    - nutation angle                 rad
%    trueeps     - true obliquity of the ecliptic rad
%    meaneps     - mean obliquity of the ecliptic rad
%    omega       -                                rad
%    nut         - matrix for tod - mod 
%    st          - matrix for pef - tod 
%    stdot       - matrix for pef - tod rate
%    pm          - matrix for ecef - pef 
%
%  coupling      :
%   nutation     - rotation for nutation          
%   sidereal     - rotation for sidereal time     
%   polarm       - rotation for polar motion      
%
%  references    :
%    vallado       2004, 219-228
%
% [rmod,vmod,amod] = ecef2mod (recef,vecef,aecef,ttt,jdut1,lod,xp,yp,eqeterms,ddpsi,ddeps )
% ----------------------------------------------------------------------------

function [rmod,vmod,amod] = ecef2mod  ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,eqeterms,ddpsi,ddeps )

        % ---- find matrices
        [deltapsi,trueeps,meaneps,omega,nut] = nutation(ttt,ddpsi,ddeps);

        [st,stdot] = sidereal(jdut1,deltapsi,meaneps,omega,lod,eqeterms );

        [pm] = polarm(xp,yp,ttt,'80');

        % ---- perform transformations
        thetasa= 7.29211514670698e-05 * (1.0  - lod/86400.0 );
        omegaearth = [0; 0; thetasa;];

%trueeps-meaneps
%deltapsi
%nut
        rpef = pm*recef;
        rmod = nut*st*rpef;

        vpef = pm*vecef;
        vmod = nut*st*(vpef + cross(omegaearth,rpef));

        temp = cross(omegaearth,rpef);
        amod = nut*st*( pm*aecef + cross(omegaearth,temp) ...
               + 2.0*cross(omegaearth,vpef) );

