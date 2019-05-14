% ----------------------------------------------------------------------------
%
%                           function setcov
%
%  this function sets the intilal states for the covariance calculations.
%
%  author        : david vallado                  719-573-2600   21 jul 2003
%
%  revisions
%
%  inputs          description                    range / units
%    reci        - position vector                km
%    veci        - velocity vector                km/s
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    terms       - number of terms for ast calculation 0,2
%
%  outputs       :
%    cartstate   - 6x1 cartesian state            m, m/s
%    classstate  - 6x1 classical orbital elements
%    flstate     - 6x1 flight orbital elements
%    eqstate     - 6x1 equinoctial orbital elements
%    fr          - retrograde factor for orbits with inlc > 90 deg  1, -1
%
%  locals        :
%
%  coupling      :
%    none
%
%  references    :
%    none
%
%  [cartstate,classstate,flstate,eqstate] = setcov(reci,veci, ...
%                                           year,mon,day,hr,min,sec,dut1,dat, ...
%                                           ttt,jdut1,lod,xp,yp,terms,printopt,anom);
% ----------------------------------------------------------------------------

function [cartstate,classstate,flstate,eqstate, fr] = setcov(reci,veci, ...
                                           year,mon,day,hr,min,sec,dut1,dat,...
                                           ttt,jdut1,lod,xp,yp,terms,printopt,anom,anom1,ddpsi,ddeps)

        constastro;
        rad = 180/pi;

        fr = 1.0;
        
        cartstate = [reci(1);reci(2);reci(3);veci(1);veci(2);veci(3)];  % in km

        % -------- convert to a classical orbit state
        [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
        classstate(1) = a*1000;
        classstate(2) = ecc;
        classstate(3) = incl;
        classstate(4) = omega;
        classstate(5) = argp;
        if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
            classstate(6) = m;
          else
            if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
                classstate(6) = nu;
            end
        end

        % -------- convert to a flight orbit state
        [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt ( reci,veci, ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );
        if strcmp(anom1,'radec') == 1  % 1 is true
            flstate(1) = rtasc; 
            flstate(2) = decl;  
          else
            if strcmp(anom1,'latlon') == 1  % 1 is true
                flstate(1) = lon;
                flstate(2) = latgc;
            end
        end
        flstate(3) = fpa;
        flstate(4) = az;
        flstate(5) = magr*1000;  % convert to m
        flstate(6) = magv*1000;

        % test position and velocity going back
        avec = [0;0;0];
        [recef,vecef,aecef] = eci2ecef(reci,veci,avec,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);
        vx = magv*( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) ); 
        vy = magv*( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );  
        vz = magv*( sin(latgc)*sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
 vtecef = [vx;  vy;  vz]'
        % correct:
        ve1 = magv*( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); % m/s
        ve2 = magv*( -sin(rtasc)*sin(decl)*cos(az)*cos(fpa) + cos(rtasc)*sin(az)*cos(fpa) + sin(rtasc)*cos(decl)*sin(fpa) );  
        ve3 = magv*( sin(decl)*sin(fpa) + cos(decl)*cos(az)*cos(fpa) );
 veci'
 vteci = [ve1;  ve2;  ve3]'

        % -------- convert to an equinoctial orbit state
        [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq ( reci, veci );
        if strcmp(anom,'meana') == 1 || strcmp(anom,'truea') == 1  % 1 is true
            eqstate(1) = a*1000;
        else
            if strcmp(anom,'meann') == 1 || strcmp(anom,'truen') == 1  % 1 is true
                eqstate(1) = n;
            end 
        end        
        eqstate(2) = af;
        eqstate(3) = ag;
        eqstate(4) = chi;
        eqstate(5) = psi;
        if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
            eqstate(6) = meanlonM;
          else
            if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
                eqstate(6) = meanlonNu;
            end
        end
        
        if printopt == 'y'
            % --------------------- write out input data --------------------------
            fprintf(1,'input data \n' );
            fprintf(1,'year %5i ',year);
            fprintf(1,'mon %4i ',mon);
            fprintf(1,'day %3i ',day);
            fprintf(1,'hr %3i:%2i:%8.6f\n',hr,min,sec );
            fprintf(1,'dut1 %8.6f s',dut1);
            fprintf(1,' dat %3i s',dat);
            fprintf(1,' xp %8.6f "',xp);
            fprintf(1,' yp %8.6f "',yp);
            fprintf(1,' lod %8.6f s\n',lod);
            fprintf(1,'r    %14.7f%14.7f%14.7f',reci );
            fprintf(1,' v %14.9f%14.9f%14.9f\n',veci );

            fprintf(1,'          p km       a km      ecc      incl deg    ');
            fprintf(1,' raan deg     argp deg      nu deg      m deg \n');
            fprintf(1,'coes %11.4f %11.4f %11.7f %11.5f %11.5f', ...
                    p,a,ecc,incl*rad,omega*rad );
            fprintf(1,'%11.5f %11.5f %11.5f\n',...
                    argp*rad,nu*rad,m*rad );

            fprintf(1,'           a           af           ag');
            fprintf(1,'           chi        psi            meanlonM     meanLonNu   fr\n');
            fprintf(1,'eq   %14.7f %14.7f %14.7f %15.7f %14.7f %14.7f %14.7f %2.0f\n',a, af, ag, chi, psi, meanlonM*rad, meanlonNu*rad, fr);

            fprintf(1,'       lon deg       latgc deg     rtasc deg      decl deg      fpa deg       ');
            fprintf(1,' az deg       magr km      magv km/s\n');
            fprintf(1,'flt  %14.7f%14.7f%14.7f%14.7f%14.7f%15.7f%14.7f%14.7f\n', ...
                    lon*rad,latgc*rad,rtasc*rad,decl*rad,fpa*rad,az*rad,magr,magv );
        end;

