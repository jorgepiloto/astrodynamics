    % ------------------------------------------------------------------------------
    %
    %                           function lambertb
    %
    %  this function solves lambert's problem using battins method. the method is
    %    developed in battin (1987) and explained by Thompson 2018. it uses continued
    %    fractions to speed the solution and has several parameters that are defined
    %    differently than the traditional gaussian technique.
    %
    %  author        : david vallado                  719-573-2600   12 feb 2018
    %
    %  inputs          description                    range / units
    %    r1          - ijk position vector 1          km
    %    v1          - ijk velocity vector 1 needed for 180 deg transfer  km / s
    %    r2          - ijk position vector 2          km
    %    dm          - dir of motion (long, short)        'l','s'
    %                  this is really a period discriminator
    %    df          - dir of flight (direct, retrograde) 'd','r'
    %                  this is the inclination discriminator
    %    nrev        - number of revs to complete     0, 1, ...
    %    dtsec       - time between r1 and r2         s
    %
    %  outputs       :
    %    v1t         - ijk velocity vector            km / s
    %    v2t         - ijk velocity vector            km / s
    %    error       - error flag                     'ok',...
    %
    %  locals        :
    %    i           - index
    %    loops       -
    %    u           -
    %    b           -
    %    sinv        -
    %    cosv        -
    %    rp          -
    %    x           -
    %    xn          -
    %    y           -
    %    l           -
    %    m           -
    %    cosdeltanu  -
    %    sindeltanu  -
    %    dnu         -
    %
    %  coupling      :
    %    mag         - magnitude of a vector
    %    arcsinh     - inverse hyperbolic sine
    %    arccosh     - inverse hyperbolic cosine
    %    sinh        - hyperbolic sine
    %
    %  references    :
    %    vallado       2013, 493-497, ex 7-5
    %    thompson      2018
    %
    % [v1dv, v2dv, errorb] = lambertb ( r1, v1, r2, dm, df, nrev, dtsec );
    % ------------------------------------------------------------------------------

    function [v1dv, v2dv, errorb] = lambertb ( r1, v1, r2, dm, df, nrev, dtsec )
    constmath;
    constastro;

    errorb = '      ok';
    y = 0.0;
    k2 = 0.0;
    u = 0.0;
    v1dv = [1000,1000,1000];
    v2dv = [1000,1000,1000];

    magr1 = mag(r1);
    magr2 = mag(r2);

    cosdeltanu= dot(r1,r2)/(magr1*magr2);
    % make sure it's not more than 1.0
    if (abs(cosdeltanu) > 1.0)
        cosdeltanu = 1.0 * sign(cosdeltanu);
    end

    rcrossr = cross( r1,r2 );
    magrcrossr = mag(rcrossr);
    if dm == 's'
        sindeltanu= magrcrossr/(magr1*magr2);
    else
        sindeltanu= -magrcrossr/(magr1*magr2);
    end;
    dnu   = atan2( sindeltanu, cosdeltanu );
    % the angle needs to be positive to work for the long way
    if dnu < 0.0
        dnu = 2.0*pi + dnu;
    end

    % these are the same
    chord= sqrt( magr1*magr1 + magr2*magr2 - 2.0*magr1*magr2*cosdeltanu );
    % chord= mag(r2-r1);

    s    = (magr1 + magr2 + chord)*0.5;
    ror   = magr2/magr1;
    eps   = ror - 1.0;

    lam = 1.0/s * sqrt(magr1*magr2) * cos(dnu*0.5);
    L = ((1.0 - lam)/(1.0 + lam))^2;
    m = 8.0*mu*dtsec*dtsec / (s^3*(1.0 + lam)^6);
    %        tan2w = 0.25*eps*eps / (sqrt(ror) + ror * (2.0 + sqrt(ror) ) );
    %        rp    = sqrt(magr1*magr2)*( (cos(dnu*0.25))^2 + tan2w );
    %        if ( dnu < pi )
    %            L = ( (sin(dnu*0.25))^2 + tan2w ) / ( (sin(dnu*0.25))^2 + tan2w + cos(dnu*0.5) );
    %        else
    %            L = ( (cos(dnu*0.25))^2 + tan2w - cos(dnu*0.5) ) / ( (cos(dnu*0.25))^2 + tan2w );
    %        end
    %        m    = mu * dtsec*dtsec / ( 8.0*rp*rp*rp );
    % initial guess
    if (nrev > 0)
        xn = 1.0 + 4.0*L;
    else
        xn   = L;   %l    % 0.0 for par and hyp, l for ell
    end

    %    lim1 = sqrt(m/L);
    % alt approach for high energy (long way, retro multi-rev) case
    if (df == 'r') && (nrev > 0)
        xn = 1e-20;  % be sure to reset this here!!
        x    = 10.0;  % starting value
        loops = 1;
        while ((abs(xn-x) >= small) && (loops <= 20))
            x = xn;
            temp = 1.0 / (2.0*(L - x*x));
            temp1 = sqrt(x);
            temp2 = (nrev*pi*0.5 + atan(temp1)) / temp1;
            h1   = temp * (L + x) * (1.0 + 2.0*x + L);
            h2   = temp * m * temp1 * ((L - x*x) * temp2 - (L + x));

            b = 0.25*27.0*h2 / ( (temp1*(1.0+h1))^3 );
            if b < -1.0 % reset the initial condition
                f = 2.0 * cos(1.0/3.0 * acos(sqrt(b + 1.0)));
            else
                A = (sqrt(b) + sqrt(b+1.0))^(1.0/3.0);
                f = A + 1.0/A;
            end

            y  = 2.0/3.0 * temp1 * (1.0 + h1) *(sqrt(b + 1.0) / f + 1.0 );
            xn = 0.5 * ((m/(y*y) - (1.0 + L)) - sqrt((m/(y*y) - (1.0 + L))^2 - 4.0*L));
            loops = loops + 1;
        end  % while
        fprintf(1,' %3i yh %11.6f x %11.6f h1 %11.6f h2 %11.6f b %11.6f f %11.7f \n',loops, y, x, h1, h2, b, f );
        x = xn;
        a = s*(1.0 + lam)^2*(1.0 + x)*(L + x) / (8.0*x);
        p = (2.0*magr1*magr2*(1.0 + x)*sin(dnu*0.5)^2) / (s*(1 + lam)^2 * (L + x));  % thompson
        ecc = sqrt(1.0 - p/a);
        [v1dv, v2dv] = lambhodograph( r1, v1, r2, p, ecc, dnu, dtsec );
        fprintf(1,'high v1t %16.8f %16.8f %16.8f \n',v1dv );
    else
        % standard processing
        % note that the dr nrev=0 case is not represented
        loops= 1;
        y1 = 0.0;
        x    = 10.0;  % starting value
        while ((abs(xn-x) >= small) && (loops <= 30))
            if (nrev > 0)
                x = xn;
                temp = 1.0 / ( (1.0 + 2.0*x + L) * (4.0*x) );
                temp1 = (nrev*pi*0.5 + atan(sqrt(x))) / sqrt(x);
                h1   = temp * (L + x)^2 * (3.0*(1.0 + x)^2 * temp1 - (3.0 + 5.0*x));
                h2   = temp * m * ((x*x - x*(1.0 + L) - 3.0*L) * temp1 + (3.0*L +x));
            else
                x    = xn;
                tempx  = seebatt(x);
                denom= 1.0 / ( (1.0 + 2.0*x + L) * (4.0*x + tempx*(3.0 + x) ) );
                h1   = (L + x)^2 * (1.0 + 3.0*x + tempx)*denom;
                h2   = m*(x - L + tempx)*denom;
            end

            % ----------------------- evaluate cubic ------------------
            b = 0.25*27.0*h2 / ((1.0+h1)^3 );

            %        if b < -1.0 % reset the initial condition
            %fprintf(1,'xx %11.6f  %11.6f  %11.6f \n',L, xn, b);
            %            xn = 1.0 - 2.0*L
            %        end
            %        else
            %            if y1 > lim1
            %                xn = xn * (lim1/y1)
            %            end
            %            else
            u  = 0.5*b / ( 1.0 + sqrt(1.0 + b) );
            k2 = kbatt(u);
            y  = ( (1.0 + h1) / 3.0 )*(2.0 + sqrt(1.0 + b) / (1.0 + 2.0*u*k2*k2) );
            xn= sqrt( ( (1.0 - L)*0.5 )^2 + m/(y*y) ) - (1.0 + L)*0.5;
            %                    xn = sqrt(l*l + m/(y*y)) - (1.0 - l); alt, doesn't seem to work
            %            end
            %        end

            y1=  sqrt( m/((L + x)*(1.0 + x)) );
            loops = loops + 1;
  %        fprintf(1,' %3i yb %11.6f x %11.6f k2 %11.6f b %11.6f u %11.6f y1 %11.7f \n',loops,y, x, k2, b, u, y1 );
        end  % while
        fprintf(1,' %3i yb %11.6f x %11.6f k2 %11.6f b %11.6f u %11.6f y1 %11.7f \n',loops,y, x, k2, b, u, y1 );

        if (loops < 30)
            % blair approach use y from solution
            %       lam = 1.0/s * sqrt(magr1*magr2) * cos(dnu*0.5);
            %       m = 8.0*mu*dtsec*dtsec / (s^3*(1.0 + lam)^6);
            %       L = ((1.0 - lam)/(1.0 + lam))^2;
            %a = s*(1.0 + lam)^2*(1.0 + x)*(lam + x) / (8.0*x);
            % p = (2.0*magr1*magr2*(1.0 + x)*sin(dnu*0.5)^2)^2 / (s*(1 + lam)^2*(lam + x));  % loechler, not right?
            p = (2.0*magr1*magr2*y*y*(1.0 + x)^2*sin(dnu*0.5)^2) / (m*s*(1 + lam)^2);  % thompson
            ecc = sqrt( (eps^2 + 4.0*magr2/magr1*sin(dnu*0.5)^2*((L-x)/(L+x))^2) / (eps^2 + 4.0*magr2/magr1*sin(dnu*0.5)^2)  );
            [v1dv, v2dv] = lambhodograph( r1, v1, r2, p, ecc, dnu, dtsec );
%            fprintf(1,'oldb v1t %16.8f %16.8f %16.8f %16.8f\n',v1dv, mag(v1dv) );
            %         r_180 = 0.001;  % 1 meter
            %         [v1dvh,v2dvh] = lambert_vel(r1, v1, r2, dnu, p, ecc, mu, dtsec, r_180);
            %         fprintf(1,'newb v1t %16.8f %16.8f %16.8f %16.8f\n',v1dvh, mag(v1dvh) );

            % Battin solution to orbital parameters (and velocities)
            % thompson 2011, loechler 1988
            if dnu > pi
                lam = -sqrt((s-chord)/s);
            else
                lam = sqrt((s-chord)/s);
            end
            %      x = xn;

            % loechler pg 21 seems correct!
            v1dvl = 1.0/(lam*(1.0 + lam))*sqrt(mu*(1.0+x)/(2.0*s^3*(L + x)))*((r2-r1) + s*(1.0+lam)^2*(L + x)/(magr1*(1.0 + x))*r1);
            % added v2
            v2dvl = 1.0/(lam*(1.0 + lam))*sqrt(mu*(1.0+x)/(2.0*s^3*(L + x)))*((r2-r1) - s*(1.0+lam)^2*(L + x)/(magr2*(1.0 + x))*r2);
            %fprintf(1,'loe v1t %16.8f %16.8f %16.8f %16.8f\n',v1dvl, mag(v1dvl) );
            %fprintf(1,'loe v2t %16.8f %16.8f %16.8f %16.8f\n',v2dvl, mag(v2dvl) );
        end  % if loops converged < 30
    end
