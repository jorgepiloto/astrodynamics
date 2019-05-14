%
% ----------------------------------------------------------------------------
%
%                           function rv2flt.m
%
%  this function transforms a position and velocity vector into the flight
%    elements - latgc, lon, fpa, az, position and velocity magnitude.
%
%  author        : david vallado                  719-573-2600   17 jun 2002
%
%  revisions
%    vallado     - add terms for ast calculation                 30 sep 2002
%    vallado     - chg magr var names                            23 may 2003
%
%  inputs          description                    range / units
%    r           - eci position vector            km
%    v           - eci velocity vector            km/s
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    terms       - number of terms for ast calculation 0,2
%    ddpsi,ddeps - corrections for fk5 to gcrf    rad
%
%  outputs       :
%    magr        - eci position vector magnitude  km
%    magv        - eci velocity vector magnitude  km/sec
%    latgc       - geocentric latitude            rad
%    lon         - longitude                      rad
%    fpa         - sat flight path angle          rad
%    az          - sat flight path az             rad
%
%  locals        :
%    fpav        - sat flight path anglefrom vert rad
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2001, xx
%
% [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ( r,v,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );
% ----------------------------------------------------------------------------

function [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ( reci,veci,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );

        twopi = 2.0*pi;

        small = 0.00000001;

        magr = mag(reci);
        magv = mag(veci);

        % -------- convert r to ecef for lat/lon calculation
        a = [0;0;0];
        [recef,vecef,aecef] = eci2ecef(reci,veci,a,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);

        % ----------------- find longitude value  ----------------- uses ecef
        temp = sqrt( recef(1)*recef(1) + recef(2)*recef(2) );
        if ( temp < small )
            lon= atan2( vecef(2), vecef(1) );
          else
            lon= atan2( recef(2), recef(1) );
        end

        %latgc = atan2( recef(3) , sqrt(recef(1)^2 + recef(2)^2) )
        latgc = asin( recef(3) / magr );

        % ------------- calculate rtasc and decl ------------------ uses eci
        temp= sqrt( reci(1)*reci(1) + reci(2)*reci(2) );
        if ( temp < small )
            rtasc= atan2( veci(2) , veci(1) );
          else
            rtasc= atan2( reci(2) , reci(1) );
        end
        %decl = atan2( reci(3) , sqrt(reci(1)^2 + reci(2)^2) )
        decl = asin( reci(3)/magr );

        h = cross(reci,veci);
        hmag = mag(h);
        rdotv= dot(reci,veci);
        fpav= atan2(hmag,rdotv);
        fpa = pi*0.5 - fpav;

        hcrossr = cross(h,reci);

        az = atan2( reci(1)*hcrossr(2) - reci(2)*hcrossr(1), hcrossr(3)*magr );

