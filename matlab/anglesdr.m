% ------------------------------------------------------------------------------
%
%                           function anglesdr
%
%  this function solves the problem of orbit determination using three
%    optical sightings.  the solution function uses the double-r technique.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%   8 oct 2007
%
%  inputs          description                    range / units
%    rtasc1       - right ascension #1            rad
%    rtasc2       - right ascension #2            rad
%    rtasc3       - right ascension #3            rad
%    decl1        - declination #1                rad
%    decl2        - declination #2                rad
%    decl3        - declination #3                rad
%    jd1          - julian date of 1st sighting   days from 4713 bc
%    jd2          - julian date of 2nd sighting   days from 4713 bc
%    jd3          - julian date of 3rd sighting   days from 4713 bc
%    rs           - ijk site position vector      km
%
%  outputs        :
%    r            - ijk position vector at t2     km
%    v            - ijk velocity vector at t2     km / s
%
%  locals         :
%    l1           - line of sight vector for 1st
%    l2           - line of sight vector for 2nd
%    l3           - line of sight vector for 3rd
%    tau          - taylor expansion series about
%                   tau ( t - to )
%    tausqr       - tau squared
%    i            - index
%    d            -
%    rho          - range from site to sat at t2  km
%    rhodot       -
%    dmat         -
%    rs1          - site vectors
%    rs2          -
%    rs3          -
%    earthrate    - velocity of earth rotation
%    p            -
%    q            -
%    oldr         -
%    oldv         -
%    f1           - f coefficient
%    g1           -
%    f3           -
%    g3           -
%    l2dotrs      -
%
%  coupling       :
%    mag          - magnitude of a vector
%    matmult      - multiply two matrices together
%    angl         - angl between two vectors
%
%  references     :
%    vallado       2007, 439-443
%
% [r2,v2] = anglesdr ( decl1,decl2,decl3,rtasc1,rtasc2,rtasc3,jd1,jd2,jd3,latgd,lon,alt );
% ------------------------------------------------------------------------------  

function [r2,v2] = anglesdr ( decl1,decl2,decl3,rtasc1,rtasc2, ...
                              rtasc3,jd1,jd2,jd3,rsite1,rsite2,rsite3, re, mu, tu );

% -------------------------  implementation   -------------------------
%   constastro;

% for sun
%re = 149597870.0;
%mu = 1.32712428e11;
   rad = 180.0 / pi;
   
   magr1in = 2.0*re;
   magr2in = 2.01*re;
   direct  = 'y';

   tol    = 1e-8*re;   % km
   pctchg = 0.005;
 
   % subtract dates and convert fraction of day to seconds
   t1 = (jd1 - jd2)*tu;  % secs
   t3 = (jd3 - jd2)*tu;

   % form line of sight vectors
   los1 = [cos(decl1)*cos(rtasc1) cos(decl1)*sin(rtasc1) sin(decl1)]';
   los2 = [cos(decl2)*cos(rtasc2) cos(decl2)*sin(rtasc2) sin(decl2)]';
   los3 = [cos(decl3)*cos(rtasc3) cos(decl3)*sin(rtasc3) sin(decl3)]';

   % --------- now we're ready to start the actual double r algorithm ---------
   magr1old  = 99999.9;
   magr2old  = 99999.9;
   magrsite1 = mag(rsite1);
   magrsite2 = mag(rsite2);
   magrsite3 = mag(rsite3);

   % take away negatives because escobal defines rs opposite
   cc1 = 2.0*dot(los1,rsite1);
   cc2 = 2.0*dot(los2,rsite2);
   ktr = 0;

   % main loop to get three values of the double-r for processing
   while (abs(magr1in-magr1old) > tol | abs(magr2in-magr2old) > tol)
       ktr = ktr + 1;
       fprintf(1,'%2i ',ktr);
       [r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
                       los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, re, mu);

       % check intermediate status
       f  = 1.0 - a/magr2*(1.0-cos(deltae32));
       g  = t3 - sqrt(a^3/mu)*(deltae32-sin(deltae32));
       v2 = (r3 - f*r2)/g;
       [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2,v2, re, mu);
       fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
              p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );

       % -------------- re-calculate f1 and f2 with r1 = r1 + delta r1
       magr1o = magr1in;
       magr1in = (1.0+pctchg)*magr1in;
       deltar1 = pctchg*magr1in;
       [r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
                              los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, re, mu);
       pf1pr1 = (f1delr1-f1)/deltar1;
       pf2pr1 = (f2delr1-f2)/deltar1;
    
       % ----------------  re-calculate f1 and f2 with r2 = r2 + delta r2
       magr1in = magr1o;
       deltar1 = pctchg*magr1in;
       magr2o = magr2in;
       magr2in = (1.0+pctchg)*magr2in;
       deltar2 = pctchg*magr2in;
       [r2,r3,f1delr2,f2delr2,q3,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
                              los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, re, mu);
       pf1pr2 = (f1delr2-f1)/deltar2;
       pf2pr2 = (f2delr2-f2)/deltar2;

%       f  = 1.0 - a/magr2*(1.0-cos(deltae32));
%       g  = t3 - sqrt(a^3/mu)*(deltae32-sin(deltae32));
%       v2 = (r3 - f*r2)/g;
%       [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2,v2);
%       fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
%              p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );


       % ------------ now calculate an update
       magr2in = magr2o;
       deltar2 = pctchg*magr2in;

       delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
       delta1 = pf2pr2*f1 - pf1pr2*f2;
       delta2 = pf1pr1*f2 - pf2pr1*f1;

       deltar1 = -delta1/delta;
       deltar2 = -delta2/delta;
     
       magr1old = magr1in;
       magr2old = magr2in;

       %  may need to limit the amount of the correction
       if abs(deltar1) > magr1in*pctchg
         fprintf(1,'%11.7f \n',deltar1 );
%         deltar1 = sign(deltar1)*magr1in*pctchg;
        end;
       if abs(deltar2) > magr2in*pctchg
         fprintf(1,'%11.7f \n',deltar2 );
%         deltar2 = sign(deltar2)*magr2in*pctchg;
        end;

       magr1in = magr1in + deltar1;
       magr2in = magr2in + deltar2;

       fprintf(1,'qs %11.7f  %11.7f  %11.7f \n',q1,q2,q3);
       fprintf(1,'magr1o %11.7f delr1 %11.7f magr1 %11.7f %11.7f  \n',magr1o,deltar1,magr1in,magr1old);
       fprintf(1,'magr2o %11.7f delr2 %11.7f magr2 %11.7f %11.7f  \n',magr2o,deltar2,magr2in,magr2old);


%       f  = 1.0 - a/magr2*(1.0-cos(deltae32));
%       g  = t3 - sqrt(a^3/mu)*(deltae32-sin(deltae32));
%       v2 = (r3 - f*r2)/g;
%       [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2,v2);
%       fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
%              p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );

       fprintf(1,'=============================================== \n');

     pause;
   end;

   % needed to get the r2 set properly since the last one was moving r2
   [r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
                          los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, re, mu);

   f  = 1.0 - a/magr2*(1.0-cos(deltae32));
   g  = t3 - sqrt(a^3/mu)*(deltae32-sin(deltae32));
   v2 = (r3 - f*r2)/g;

