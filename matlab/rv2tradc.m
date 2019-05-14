%
% ------------------------------------------------------------------------------
%
%                           function rv2tradc
%
%  this function converts geocentric equatorial (eci) position and velocity
%    vectors into range, topcentric right acension, declination, and rates.  
%    notice the value of small as it can affect the rate term calculations. 
%    the solution uses the velocity vector to find the singular cases. also,
%    the right acension and declination rate terms are not observable unless 
%    the acceleration vector is available.
%
%  author        : david vallado                  719-573-2600   19 jul 2004
%
%  revisions
%
%  inputs          description                    range / units
%    reci        - eci position vector            km
%    veci        - eci velocity vector            km/s
%    latgd       - geodetic latitude              -pi/2 to pi/2 rad
%    lon         - longitude of site              -2pi to 2pi rad
%    alt         - altitude                       km
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       rad
%    yp          - polar motion coefficient       rad
%    terms       - number of terms for ast calculation 0,2
%
%  outputs       :
%    rho         - satellite range from site      km
%    rtasc       - topocentric right ascension    0.0 to 2pi rad
%    decl        - topocentric declination        -pi/2 to pi/2 rad
%    drho        - range rate                     km/s
%    daz         - xxazimuth rate                   rad / s
%    del         - xxelevation rate                 rad / s
%
%  locals        :
%    rhoveci     - eci range vector from site     km
%    drhoveci    - eci velocity vector from site  km / s
%    rhoeci      - eci range vector from site     km
%    drhoeci     - sez velocity vector from site  km
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
%    vallado       2001, 250-255, alg 27
%
% [rho,trtasc,tdecl,drho,dtrtasc,dtdecl] = rv2tradc ( reci,veci, latgd,lon,alt,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );
% ------------------------------------------------------------------------------

function [rho,trtasc,tdecl,drho,dtrtasc,dtdecl] = rv2tradc ( reci,veci, latgd,lon,alt,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );

        constmath;

        % --------------------- implementation ------------------------
        % ----------------- get site vector in ecef -------------------
        [rs,vs] = site ( latgd,lon,alt );

%rs
%vs
        % -------------------- convert ecef to eci --------------------
        a = [0;0;0];
        [rseci,vseci,aeci] = ecef2eci(rs,vs,a,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps);
%rseci
%vseci

%rseci = rs;
%vseci = vs;
%[recef,vecef,aecef] = eci2ecef(reci,veci,aeci,ttt,jdut1,lod,xp,yp,2,0,0);
%reci = recef;
%veci = vecef;

        % ------- find eci range vector from site to satellite -------
        rhoeci  = reci - rseci;
        drhoeci = veci - vseci;
        rho      = mag(rhoeci);

        % ------------- calculate azimuth and elevation ---------------
        temp= sqrt( rhoeci(1)*rhoeci(1) + rhoeci(2)*rhoeci(2) );
        if ( temp < small )
            trtasc = atan2( drhoeci(2), drhoeci(1) );
          else
            trtasc= atan2( rhoeci(2), rhoeci(1) );
          end

        if ( ( temp < small ) )           % directly over the north pole
            tdecl= sign(rhoeci(3))*halfpi;   % +- 90 deg
          else
            magrhoeci = mag(rhoeci);
            tdecl= asin( rhoeci(3) / magrhoeci );
          end

        % ------ calculate range, azimuth and elevation rates ---------
        temp1= -rhoeci(2)*rhoeci(2) - rhoeci(1)*rhoeci(1);
        drho= dot(rhoeci,drhoeci)/rho;
        if ( abs( temp1 ) > small )
            dtrtasc= ( drhoeci(1)*rhoeci(2) - drhoeci(2)*rhoeci(1) ) / temp1;
          else
            dtrtasc= 0.0;
          end

        if ( abs( temp ) > small )
            dtdecl= ( drhoeci(3) - drho*sin( tdecl ) ) / temp;
          else
            dtdecl= 0.0;
          end

