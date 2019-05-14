%
% ------------------------------------------------------------------------------
%
%                           function pkepler
%
%  this function propagates a satellite's position and velocity vector over
%    a given time period accounting for perturbations caused by j2.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    ro          - original position vector       km
%    vo          - original velocity vector       km/sec
%    ndot        - time rate of change of n       rad/sec
%    nddot       - time accel of change of n      rad/sec2
%    dtsec       - change in time                 sec
%
%  outputs       :
%    r           - updated position vector        km
%    v           - updated velocity vector        km/sec
%
%  locals        :
%    p           - semi-paramter                  km
%    a           - semior axis                    km
%    ecc         - eccentricity
%    incl        - inclination                    rad
%    argp        - argument of periapsis          rad
%    argpdot     - change in argument of perigee  rad/sec
%    omega       - longitude of the asc node      rad
%    omegadot    - change in omega                rad
%    e0          - eccentric anomaly              rad
%    e1          - eccentric anomaly              rad
%    m           - mean anomaly                   rad/sec
%    mdot        - change in mean anomaly         rad/sec
%    arglat      - argument of latitude           rad
%    arglatdot   - change in argument of latitude rad/sec
%    truelon     - true longitude of vehicle      rad
%    truelondot  - change in the true longitude   rad/sec
%    lonper     - longitude of periapsis         rad
%    lonperodot  - longitude of periapsis change  rad/sec
%    n           - mean angular motion            rad/sec
%    nuo         - true anomaly                   rad
%    j2op2       - j2 over p sqyared
%    sinv,cosv   - sine and cosine of nu
%
%  coupling:
%    rv2coe      - orbit elements from position and velocity vectors
%    coe2rv      - position and velocity vectors from orbit elements
%    newtonm     - newton rhapson to find nu and eccentric anomaly
%
%  references    :
%    vallado       2007, 687, alg 64
%
% [r,v] = pkepler( ro,vo, dtsec, ndot,nddot );
% ----------------------------------------------------------------------------- }

function [r,v] = pkepler( ro,vo, dtsec, ndot,nddot );

     constastro;
     j2 =  0.00108263;

     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo);
%     fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
%     fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n',...
%             p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
%             arglat*rad,truelon*rad,lonper*rad );

     n= sqrt(mu/(a*a*a));

     % ------------- find the value of j2 perturbations -------------
     j2op2   = (n*1.5*re^2*j2) / (p*p);
%     nbar    = n*( 1.0 + j2op2*sqrt(1.0-ecc*ecc)* (1.0 - 1.5*sin(incl)*sin(incl)) );
     omegadot= -j2op2 * cos(incl);
     argpdot =  j2op2 * (2.0-2.5*sin(incl)*sin(incl));
     mdot    =  n;

     a     = a - 2.0*ndot*dtsec * a / (3.0*n);
     ecc   = ecc - 2.0*(1.0 - ecc)*ndot*dtsec / (3.0*n);
     p     = a*(1.0 - ecc*ecc);

     % ----- update the orbital elements for each orbit type --------
     if ecc < small
         % -------------  circular equatorial  ----------------
         if ( incl < small ) || ( abs(incl-pi) < small )
             truelondot= omegadot + argpdot + mdot;
             truelon   = truelon  + truelondot * dtsec;
             truelon   = rem(truelon, twopi);
           else
           % -------------  circular inclined    --------------
              omega    = omega + omegadot * dtsec;
              omega    = rem(omega, twopi);
              arglatdot= argpdot + mdot;
              arglat   = arglat + arglatdot * dtsec;
              arglat   = rem(arglat, twopi);
        end;
        else
          % -- elliptical, parabolic, hyperbolic equatorial ---
          if ( incl < small ) || ( abs(incl-pi) < small )
              lonperdot= omegadot + argpdot;
              lonper  = lonper + lonperdot * dtsec;
              lonper  = rem(lonper, twopi);
              m        = m + mdot*dtsec + ndot*dtsec*dtsec + nddot*dtsec*dtsec*dtsec;
              m        = rem(m, twopi);
              [e0,nu]= newtonm(ecc,m);
            else
            % --- elliptical, parabolic, hyperbolic inclined --
              omega= omega + omegadot * dtsec;
              omega= rem(omega, twopi);
              argp = argp  + argpdot  * dtsec;
              argp = rem(argp, twopi);
              m    = m + mdot*dtsec + ndot*dtsec*dtsec + nddot*dtsec*dtsec*dtsec;
              m    = mod(m, twopi);
              [e0,nu]= newtonm(ecc,m);
         end;
     end;

        % ------------- use coe2rv to find new vectors ---------------
        [r,v] = coe2rv(p,ecc,incl,omega,argp,nu,arglat,truelon,lonper);
        r = r';
        v = v';
%        fprintf(1,'r    %15.9f%15.9f%15.9f',r );
%        fprintf(1,' v %15.10f%15.10f%15.10f\n',v );

