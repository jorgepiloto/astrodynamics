% ------------------------------------------------------------------------------
%
%                           function findtof
%
%  this function finds the time of flight given the initial position vectors,
%    semi-parameter, and the sine and cosine values for the change in true
%    anomaly.  the result uses p-iteration theory to analytically find the result.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%    vallado     - fix tolerances                                 5 sep 2002
%
%  inputs          description                    range / units
%    ro          - interceptor position vector    km
%    r           - target position vector         km
%    p           - semiparameter                  km
%
%  outputs       :
%    tof         - time for transfer              sec
%
%  locals        :
%    sindnu      - sine of change in nu           rad
%    cosdnu      - cosine of change in nu         rad
%    deltae      -
%    deltah      -
%    k           -
%    l           -
%    m           -
%    a           -
%    f           -
%    g           -
%    fdot        -
%    sindeltae   - sine value
%    cosdeltae   - cosine value
%    rcrossr     - cross product of two positions
%
%  coupling      :
%    cross       - cross product of two vectors
%    sinh        - hyperbolic sine
%    arccosh     - arc hyperbolic cosine
%
%  references    :
%    vallado       2007, 134-135, alg 11
%
% [tof] = findtof ( ro,r, p );
% ------------------------------------------------------------------------------

function [tof] = findtof ( ro,r, p );

        small = 0.00000001;
        constastro;
        
        % -------------------------  implementation   -----------------
        magr = mag(r);
        magro = mag(ro);
        cosdnu= dot(ro,r)/(magro*magr);
        rcrossr = cross( ro,r );
        sindnu= mag(rcrossr)/(magro*magr);

        k= magro * magr*( 1.0 -cosdnu );
        l= magro + magr;
        m= magro * magr*( 1.0 +cosdnu );
        a= (m*k*p) / ((2.0 *m-l*l)*p*p + 2.0 *k*l*p - k*k);

        % ------  use f and g series to find velocity vectors  --------
        f = 1.0  - ( magr/p )*(1.0 -cosdnu);
        g = magro*magr*sindnu/sqrt(mu * p);
        alpha= 1.0 /a;

        if ( alpha > small )
            % ------------------------ elliptical ---------------------
            dnu  = atan2( sindnu,cosdnu );
            fdot = sqrt(mu /p) * tan(dnu*0.5 )* ...
                    ( ((1.0 -cosdnu)/p)-(1.0 /magro)-(1.0 /magr) );
            cosdeltae= 1.0 -(magro/a)*(1.0 -f);
            sindeltae= (-magro*magr*fdot)/sqrt(mu * a);
            deltae   = atan2( sindeltae,cosdeltae );
            tof      = g + sqrt(a*a*a/mu)*(deltae-sindeltae);
          else
            % ------------------------ hyperbolic ---------------------
            if ( alpha < -small )
                deltah = arccosh( 1.0 -(magro/a)*(1.0 -f) );
                tof    = g + sqrt(-a*a*a/mu)*(sinh(deltah)-deltah);
              else
                % -------------------- parabolic ----------------------
                dnu= atan2( sindnu,cosdnu );
                c  = sqrt( magr*magr+magro*magro - ...
                           2.0 *magr*magro*cos(dnu) );
                s  = (magro+magr+c ) * 0.5;
                tof= ( 2.0 /3.0  ) * sqrt(s*s*s*0.5/mu ) * ...
                           (1.0  -  ((s-c)/s)^1.5  );
              end
          end

