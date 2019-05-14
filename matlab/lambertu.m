% ------------------------------------------------------------------------------
%
%                           function lambertu
%
%  this function solves the lambert problem for orbit determination and returns
%    the velocity vectors at each of two given position vectors.  the solution
%    uses universal variables for calculation and a bissection technique
%    updating psi.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    r1          - ijk position vector 1          km
%    v1          - ijk velocity vector 1 needed for 180 deg transfer  km / s
%    r2          - ijk position vector 2          km
%    dm          - direction of motion            'l','s'  long short period
%    df          - direction of flight            'd','r'  direct retrograde
%    dtsec       - time between r1 and r2         s
%    nrev        - multiple revoluions            0, 1, ...
%    tbi         - time of the bottom interval - only needed for multi-rev cases
%                  this is a two-dimension array of psi and tof
%  outputs       :
%    v1          - ijk velocity vector            km / s
%    v2          - ijk velocity vector            km / s
%    error       - error flag                     'ok', ...
%
%  locals        :
%    vara        - variable of the iteration,
%                  not the semi-axis
%    y           - area between position vectors
%    upper       - upper bound for z
%    lower       - lower bound for z
%    cosdeltanu  - cosine of true anomaly change  rad
%    f           - f expression
%    g           - g expression
%    gdot        - g dot expression
%    x           - old universal variable x
%    xcubed      - x cubed
%    zold        - old value of z
%    znew        - new value of z
%    c2          - c2(z) function
%    c3          - c3(z) function
%    timenew     - new time                       s
%    small       - tolerance for roundoff errors
%    i, j        - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    dot         - dot product of two vectors
%    findc2c3    - find c2 and c3 functions
%
%  references    :
%    vallado       2013, 489-493, alg 58, ex 7-5
%
% [vo,v,errorl] = lambertu ( r1, v1, r2, dm, df, nrev, dtsec, tbi, outfile );
% ------------------------------------------------------------------------------

%function [v1dv,v2dv,errorl] = lambertu ( r1, v1, r2, dm, df, nrev, dtsec );
function [v1dv,v2dv,errorl] = lambertu ( r1, v1, r2, dm, df, nrev, dtwait, dtsec, tbi, outfile )

% -------------------------  implementation   -------------------------
        mu = 398600.4418;  % km/s^2
        small = 0.00001; % can affect cases where znew is multiples of 2pi^2
        numiter= 20;
        errorl  = '      ok';
        for i= 1 : 3
            v1dv(i) = 0.0;
            v2dv(i) = 0.0;
        end
        
        % try canonical units for testing
%constastro;
%mu = 1.0;
%r1 = r1/re;
%r2 = r2/re;
%dtsec = dtsec / tusec;

        % ---- find parameters that are constant for the initial geometry
        magr1 = mag(r1);
        magr2 = mag(r2);

        % this value stays constant in all calcs, vara changes with df
        cosdeltanu = dot(r1,r2) / (magr1*magr2);  
        if abs(cosdeltanu) > 1.0
            cosdeltanu = sign(cosdeltanu) * 1.0;
        end
        if ( df == 'r' )  %dm == 'l'
            vara = -sqrt( magr1*magr2*(1.0 + cosdeltanu) );
        else
            vara =  sqrt( magr1*magr2*(1.0 + cosdeltanu) );
        end
        
        % setup variables for speed
        oomu = 1.0 / sqrt(mu);

        % --------- set up initial bounds for the bissection ----------
        if ( nrev == 0 )  
            lower = -4.0*pi*pi;  % allow hyperbolic and parabolic solutions
            upper =  4.0*pi*pi;  % could be negative infinity for all cases
        else
            % set absolute limits for multi-rev cases
            lower = 4.0*nrev^2*pi*pi;
            upper = 4.0*(nrev + 1.0)^2*pi*pi;       
            % adjust based on long or short way if dm == 'l'
            %if ((dm == 'l') && (df == 'd')) || ((dm == 's') && (df == 'r'))
            %if ((df == 'r') && (dm == 's')) || ((df == 'd') && (dm == 'l'))      
            if ((df == 'r') && (dm == 'l')) || ((df == 'd') && (dm == 'l'))      
              upper = tbi(nrev,1);
            else
              lower = tbi(nrev,1);
            end
        end

        % ---------------  form initial guesses   ---------------------
        dtdpsi = 0.0;
        x      = 0.0;
        psinew = 0.0;
        if (nrev == 0)
            % use log to get initial guess
            % empirical relation here from 10000 random draws
            % 10000 cases up to 85000 dtsec  0.11604050x + 9.69546575
            psiold = (log(dtsec) - 9.61202327)/0.10918231;
            if psiold > upper
                psiold = upper - pi;
            end
        else
            if (df == 'd')  % dm == 's'
                psiold = lower + (upper - lower)*0.3;
            else
                psiold = lower + (upper - lower)*0.6;
            end
        end
        
        [c2,c3] = findc2c3( psiold );

        % -------  determine if  the orbit is possible at all ---------
        if ( abs( vara ) > 0.2 )   % not exactly zero
            loops  = 0;
            ynegktr= 1;  % y neg ktr
            dtnew = -10.0;
            while ((abs(dtnew-dtsec) >= small) && (loops < numiter) && (ynegktr <= 10))
                % fprintf(1,'%3i  dtnew-dtsec %11.7f yneg %3i \n',loops,dtnew-dtsec,ynegktr );
                if ( abs(c2) > small )
                    y= magr1 + magr2 - ( vara*(1.0 - psiold*c3)/sqrt(c2) );
                else
                    y= magr1 + magr2;
                end
                % ----------- check for negative values of y ----------
                if ( (vara > 0.0) && (y < 0.0) )  % ( vara > 0.0 ) &
                    ynegktr= 1;
                    while (( y < 0.0 ) && ( ynegktr < 10 ))
                        psinew = 0.8*(1.0 / c3)*( 1.0 - (magr1 + magr2)*sqrt(c2)/vara  );  
                        % -------- find c2 and c3 functions -----------
                        [c2,c3] = findc2c3( psinew );
                        psiold = psinew;
                        lower  = psiold;
                        if ( abs(c2) > small )
                            y= magr1 + magr2 - ( vara*(1.0-psiold*c3) / sqrt(c2) );
                        else
                            y= magr1 + magr2;
                        end
                        % outfile
                        fprintf(1,'yneg %3i  y %11.7f lower %11.7f c2 %11.7f psinew %11.7f yneg %3i \n',loops,y,lower,c2,psinew,ynegktr );
                        ynegktr = ynegktr + 1;
                    end % while
                end  % if  y neg

                if ( ynegktr < 10 )  
                    if ( abs(c2) > small )  
                        x= sqrt( y / c2 );
                    else
                        x= 0.0;
                    end
                    xcubed= x^3;
                   
                    dtnew    = (xcubed*c3 + vara*sqrt(y)) * oomu;
                    % try newton rhapson iteration to update psi
                    if abs(psiold) > 1e-5 
                        c2dot = 0.5/psiold * (1.0 - psiold*c3 - 2.0*c2);
                        c3dot = 0.5/psiold * (c2 - 3.0*c3);
                    else  % case for parabolic orbit
                        c2dot = -1.0/factorial(4) + 2.0*psiold/factorial(6) - 3.0*psiold^2/factorial(8) + 4.0*psiold^3/factorial(10) - 5.0*psiold^4/factorial(12);
                        c3dot = -1.0/factorial(5) + 2.0*psiold/factorial(7) - 3.0*psiold^2/factorial(9) + 4.0*psiold^3/factorial(11) - 5.0*psiold^4/factorial(13);
                    end    
                    dtdpsi = (xcubed*(c3dot - 3.0*c3*c2dot/(2.0*c2)) + 0.125*vara * (3.0*c3*sqrt(y)/c2 + vara/x)) * oomu;
                    % Newton iteration test to see if it keeps within the bounds
                    psinew =  psiold - (dtnew - dtsec)/dtdpsi;  
                    
                    % check if newton guess for psi is outside bounds (too steep a slope)
                    if abs(psinew) > upper || psinew < lower 
                        % --------  readjust upper and lower bounds -------
                        if ( dtnew < dtsec )
                            if psiold > lower
                                lower= psiold; 
                            end
                        end
                        if ( dtnew > dtsec )
                            if psiold < upper
                                upper= psiold;
                            end
                        end
                        psinew= (upper+lower) * 0.5;
                    end

                    % ------------- find c2 and c3 functions ----------
                    [c2,c3] = findc2c3( psinew );
 if nrev > 0
                     fprintf(1,'%3i  y %11.7f x %11.7f %11.7f dtnew %11.7f %11.7f %11.7f psinew %11.7f %11.7f \n', ...
                             loops,y,x,dtsec, dtnew, lower, upper, psinew, dtdpsi); %(dtnew - dtsec)/dtdpsi );  % c2dot, c3dot
 end
                    psiold = psinew;
                    loops = loops + 1;

                    % --- make sure the first guess isn't too close ---
                    if ( (abs(dtnew - dtsec) < small) && (loops == 1) );
                        dtnew= dtsec - 1.0;
                    end
                end  % if  ynegktr < 10
                
%              fprintf(1,'#%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y,x,dtnew,psinew );
%              fprintf(1,'%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,x/sqrt(re),dtnew/tusec,psinew );
%              fprintf(1,'%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,x/sqrt(re),dtnew/60.0,psinew );
            end % while loop

            if ( (loops >= numiter) || (ynegktr >= 10) )
                errorl= strcat('gnotconv',num2str(abs(dtnew - dtsec)));
                if ( ynegktr >= 10 )
                    errorl= 'y negati';
                end
            else
                % --- use f and g series to find velocity vectors -----
                f   = 1.0 - y/magr1;
                gdot= 1.0 - y/magr2;
                g   = 1.0 / (vara*sqrt( y/mu ));  % 1 over g
%fprintf(1,'%11.7f  %11.7f  %11.7f \n',f, gdot, g);
                
           %  fdot = sqrt(mu*y)*(-magr2-magr1 + y)/(magr1*magr2*vara);
           %  f*gdot - fdot*g                
                for i= 1 : 3
                    v1dv(i) = ( r2(i) - f*r1(i) )*g;
                    v2dv(i) = ( gdot*r2(i) - r1(i) )*g;
                end
            end   % if  the answer has converged
        else
            errorl= 'impos180';
 
              % do hohmann but in 3-d...
              % can't do bissection because w series is not accurate
              mum = 3.986004418e5;   % 14 m3/s2
              atx = (mum*(dtsec/(1.0*pi))^2)^(1.0/3.0);  % 1pi since half period
              v1tmag = sqrt(2.0*mum/(magr1) - mum/atx);
              v2tmag = sqrt(2.0*mum/(magr2) - mum/atx);
              wx = cross(r1,v1);
              wxu = unit(wx);
              v1dir = cross(r1,wxu); % get retro direc
              v2dir = cross(r2,wxu); % get retro direc
              v1diru = unit(v1dir);
              v2diru = unit(v2dir);
              v1t = -v1tmag*v1diru;
              v2t = -v2tmag*v2diru;
              fprintf(1,'%11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',r1, v1t);
              fprintf(1,'%11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',r2, v2t);
              
              v1dv = v1t;% / velkmps;
              v2dv = v2t;% / velkmps;

              tof = dtsec;
            pause;
%             % use JGCD 2011 v34 n6 1925 to solve 180 deg case
%             p = 2.0*magr1*magr2 / (magr1 + magr2); 
%             ecc = sqrt(1.0 - 4.0*magr1*magr2 / ((magr1+magr2)^2) );
%             dt = sqrt(pi^2/mu * (p / (1.0 - ecc^2))^3);
%             dnu = acos(cosdeltanu);
%      fprintf(1,'hodo %11.6f   %11.6f  %11.6f  %11.6f %14.10f \n',p, ecc, dt, dtsec, vara);
%             if abs(dt-dtsec) < 160
%                 [v1dv, v2dv] = lambhodograph( r1, v1, r2, p, ecc, dnu, dtsec );
%             end;    
        end  % if  var a > 0.0

        if (strcmp(errorl, '      ok') ~= 0)
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r1,v1dv);
            fprintf(outfile,'%10s %3i %3i %2s %2s %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f case %11.7f %11.7f %11.7f %11.7f %11.7f ', ...
                     errorl, loops, nrev, dm, df, dtwait, dtnew, y, x, v1dv(1), v1dv(2), v1dv(3), v2dv(1), v2dv(2), v2dv(3), lower, upper, psinew, dtdpsi, ecc); %(dtnew - dtsec)/dtdpsi, ecc );  % c2dot, c3dot
            fprintf(1,'C%3i %3i %2s %2s %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f  %11.7f \n', ...
                   loops, nrev, dm, df, dtnew, magr1, magr2, vara, y, x, psinew)
        else
            fprintf(outfile,'#%s \n',errorl);
            fprintf(1,'#%s \n',errorl);
        end
   
end         
       
        
        
       
       
