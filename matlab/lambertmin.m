%    test case
%    ro = [ 2.500000    0.000000    0.000000]*6378.137;
%    r  = [ 1.9151111   1.6069690   0.000000]*6378.137;
%
% lambert minimum energy, not min time!!
% min time is approximated by the parabolic case which is sort of a limit of
% what could be done
%
%  inputs          description                    range / units
%    r1          - ijk position vector 1          km
%    r2          - ijk position vector 2          km
%    dm          - direction of motion            'l','s'
%    nrev        - multiple revoluions            0, 1, ...
%
%  outputs       :
%    v           - ijk velocity vector            km / s
%    aminenergy  - semimajor axis min energy      km / s
%    tminenergy  - time min energy                km / s
%    tminabs     - time min energy - parabolic    km / s
%
%
function [v, aminenergy, tminenergy, tminabs] = lambertmin ( r1, r2, df, nrev )
        mu = 398600.4418;  % gravitational parameter km/s^2
         
        % ---- find parameters that are constant for the initial geometry
        magr1 = mag(r1);
        magr2 = mag(r2);
        cosdeltanu = dot(r1,r2) / (magr1*magr2);
             
        c = sqrt(magr1^2 + magr2^2 - 2.0*magr1*magr2*cosdeltanu);
        s = 0.5 * (magr1 + magr2 + c);
        aminenergy = 0.5*s;

        alphae = pi;
        betae = 2.0*asin( sqrt((s-c)/s) );
        
        if (df == 'd')
            tminenergy = sqrt(aminenergy^3/mu)*(2.0*nrev*pi + alphae - (betae-sin(betae)));
            %tminabs1 = sqrt(s^3/(8.0*mu))*(alphae - (betae-sin(betae)))
        else
            tminenergy = sqrt(aminenergy^3/mu)*(2.0*nrev*pi + alphae + (betae-sin(betae)));
            %tminabs1 = sqrt(s^3/(8.0*mu))*(alphae + (betae-sin(betae)))
        end
        % find parabolic tof - this will be the minimum limit for tof
        % negative sign should be smallest
        tminabs = 1.0/3.0*sqrt(2.0/mu)*(s^1.5-(s-c)^1.5);
        %tminabs = 1.0/3.0*sqrt(2.0/mu)*(s^1.5+(s-c)^1.5);

        % if calc min velocity
        rcrossr = cross( r1,r2 );
        magrcrossr = mag(rcrossr);
        pmin = magr1*magr2/c*(1.0 - cosdeltanu);
        if df == 'd'
            sindeltanu= magrcrossr/(magr1*magr2);
        else
            sindeltanu= -magrcrossr/(magr1*magr2);            
        end;
        for i= 1 : 3
            v(i) = sqrt(mu*pmin)/(magr1*magr2*sindeltanu)*(r2(i) - (1.0-magr2/pmin*(1.0-cosdeltanu))*r1(i));
        end
end         


