% -------------------------------------------------------------------------
% find the minimum psi values for the universal variable lambert problem
% for multi-rev cases
%  inputs          description                    range / units
%    r1          - ijk position vector 1             km
%    r2          - ijk position vector 2             km
%    dm          - direction of motion (long/short)  'l','s'
%    you could maybe pass in the max number of revs you're willing to
%    consider...
%    nrev        - multiple revoluions                0, 1, ...
%
%  outputs       :
%    tbi         - 2D array containing the dt (sec), psi (unitless) values for the minimums of
%    nrev = 1, 2, 3, ...
%

function [psib, tof] = lambertumins( r1, r2, nrev, df )
    small = 0.00000001;
    mu = 398600.4418;  % km/s^2
    oomu = 1.0 / sqrt(mu);  % for speed
    sqrtmu = sqrt(mu);
    numiter = 10; % arbitrary limit here - doens't seem to break it. 
    
    % ---- find parameters that are constant for the intiial geometry
    magr1 = mag(r1);
    magr2 = mag(r2);

    cosdeltanu = dot(r1,r2)/(magr1*magr2);  
    if ( df == 'r' )
        vara = -sqrt( magr1*magr2*(1.0+cosdeltanu) );
    else
        vara =  sqrt( magr1*magr2*(1.0+cosdeltanu) );
    end

    % ------------ find the minimum time for a nrev transfer --------------
 %   nrev = 0;
 %   for zz = 0: 4
 %       nrev = nrev + 1;
        % ---- get outer bounds for each nrev case
        lower = 4.0*nrev^2*pi*pi;
        upper = 4.0*(nrev + 1.0)^2*pi*pi;

        % ---- streamline since we know it's near the center
        upper = lower + (upper - lower)*0.6;             
        lower = lower + (upper - lower)*0.3;             

        % ---- initial estimate, just put in center 
        psiold = (upper + lower) * 0.5;
        [c2,c3] = findc2c3( psiold );

        loops = 0;
        dtdpsi = 200.0;
        while ((abs(dtdpsi) >= 0.1) && (loops < numiter) )
            if ( abs(c2) > small )
                y = magr1 + magr2 - ( vara*(1.0 - psiold*c3)/sqrt(c2) );
            else
                y = magr1 + magr2;
            end
            if ( abs(c2) > small )
                x = sqrt( y / c2 );
            else
                x = 0.0;
            end
            sqrty = sqrt(y);
            if abs(psiold) > 1e-5
                c2dot = 0.5/psiold * (1.0 - psiold*c3 - 2.0*c2);
                c3dot = 0.5/psiold * (c2 - 3.0*c3);
                c2ddot = 1.0/(4.0*psiold^2) * ((8.0-psiold)*c2 + 5.0*psiold*c3 - 4.0);
                c3ddot = 1.0/(4.0*psiold^2) * ((15.0-psiold)*c3 - 7.0*c2 + 1.0);
            else
                c2dot = -2.0/factorial(4) + 2.0*psiold/factorial(6) - 3.0*psiold^2/factorial(8) + 4.0*psiold^3/factorial(10) - 5.0*psiold^4/factorial(12);
                c3dot = -1.0/factorial(5) + 2.0*psiold/factorial(7) - 3.0*psiold^2/factorial(9) + 4.0*psiold^3/factorial(11) - 5.0*psiold^4/factorial(13);
                c2ddot = 0.0;
                c3ddot = 0.0;
            end
            % now solve this for dt = 0.0
            dtdpsi = x^3*(c3dot - 3.0*c3*c2dot/(2.0*c2))* oomu + 0.125*vara/sqrt(mu) * (3.0*c3*sqrty/c2 + vara/x);

            q = 0.25*vara*sqrt(c2) - x^2*c2dot;
            s1 = -24.0 * q * x^3 * c2 * sqrty * c3dot;
            s2 = 36.0 * q * x^3 * sqrty * c3 * c2dot - 16.0 * x^5 * sqrty * c3ddot * c2^2;
            s3 = 24.0 * x^5 * sqrty * (c3dot * c2dot*c2 + c3 * c2ddot * c2 - c3*c2dot^2) - 6.0 * vara * c3dot * y * c2 * x^2;
            s4 = -0.75 * vara^2*c3*c2^1.5*x^2 + 6.0 * vara * c3 *y*c2dot*x^2 + ( vara^2 * c2*(0.25*vara*sqrt(c2) - x^2*c2))*sqrty / x; % C(z)??

            dtdpsi2 = -(s1 + s2 + s3 + s4)/(16.0*sqrtmu*(c2^2*sqrty*x^2));
            % NR update
            psinew = psiold - dtdpsi / dtdpsi2;  

 %           fprintf(1,' %3i %12.4f %12.4f %12.4f %12.4f %11.4f %12.4f %12.4f %11.4f %11.4f \n',loops, y, dtnew, psiold, psinew, psinew - psiold, dtdpsi, dtdpsi2, lower, upper );
            psiold = psinew;
            [c2,c3] = findc2c3( psiold );
            % don't need for the loop iterations
            % dtnew = (x^3*c3 + vara*sqrty) * oomu;
            loops = loops + 1;
        end  % while
        % calculate once at the end
        dtnew = (x^3*c3 + vara*sqrty) * oomu;
        tof = dtnew;
        psib = psinew;
 %       fprintf(1,' nrev %3i  dtnew %12.5f psi %12.5f  lower %10.3f upper %10.3f %10.6f %10.6f \n',nrev, dtnew, psiold, lower, upper, c2, c3);
 %   end % for checking multi rev cases

end

