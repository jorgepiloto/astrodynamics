% ----------------------------------------------------------------------------
%
%                           function covcl2eq
%
%  this function transforms a six by six covariance matrix expressed in
%    classical elements into one expressed in equinoctial elements.
%
%  author        : david vallado                  719-573-2600   14 jul 2002
%
%  revisions
%    vallado     - major update                                  26 aug 2015
%
%  inputs          description                    range / units
%    classcov    - 6x6 classical covariance matrix
%    classstate  - 6x1 classical orbit state      (a e i O w nu/m)
%    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
%    fr          - retrograde factor               +1, -1
%
%  outputs       :
%    eqcov       - 6x6 equinoctial covariance matrix
%    tm          - transformation matrix
%
%  locals        :
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omaga       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    m           - mean anomaly                   0.0  to 2pi rad
%    mum         - gravitational paramater        m^3/s^2 NOTE Meters!
%
%  coupling      :
%    constastro
%
%  references    :
%    Vallado and Alfano 2015
%
% [eqcov, tm] = covcl2eq ( classcov, classstate, anom, fr);
% ----------------------------------------------------------------------------
 
function [eqcov, tm] = covcl2eq ( classcov, classstate, anom, fr)

        % -------- define the gravitational constant
        constastro;

        % --------- determine which set of variables is in use ---------
        % -------- parse the orbit state
        a = classstate(1);  % in m
        n = sqrt(mum/a^3);
        ecc     = classstate(2);
        incl    = classstate(3);
        omega   = classstate(4);
        argp    = classstate(5);
        if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
            m = classstate(6);
          else
            if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
                nu = classstate(6);
            end
        end
  
        % ---- partials of a wrt (a e i O w nu/m/tau)
        if strcmp(anom,'truea') == 1 || strcmp(anom,'meana') == 1  % 1 is true
            tm(1,1) =  1.0;
        else
            if strcmp(anom,'truen') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                %tm(1,1) = 1.0;
                tm(1,1) =  -( 3.0*sqrt(mum/a^3) ) / (2.0*a);  % if class = a, equin = n
            end
        end    
        
        % ---------------- calculate matrix elements ------------------
        tm(1,2) =  0.0;
        tm(1,3) =  0.0;
        tm(1,4) =  0.0;
        tm(1,5) =  0.0;
        tm(1,6) =  0.0;

        % ---- partials of af wrt (a e i O w nu/m/tau)
        tm(2,1) = 0.0;
        tm(2,2) = cos(fr*omega + argp);
        tm(2,3) = 0.0;
        tm(2,4) = -ecc*fr*sin(fr*omega + argp);
        tm(2,5) = -ecc*sin(fr*omega + argp);
        tm(2,6) = 0.0;

        % ---- partials of ag wrt (a e i O w nu/m/tau)
        tm(3,1) =  0.0;
        tm(3,2) =  sin(fr*omega + argp);
        tm(3,3) =  0.0;
        tm(3,4) =  ecc*fr*cos(fr*omega + argp);
        tm(3,5) =  ecc*cos(fr*omega + argp);
        tm(3,6) =  0.0;

        % ---- partials of chi wrt (a e i O w nu/m/tau)
        tm(4,1) =  0.0;
        tm(4,2) =  0.0;
        tm(4,3) = sin(omega) * (0.5*tan(incl*0.5)*tan(incl*0.5) + 0.5) * fr * tan(incl*0.5)^(fr-1);
        tm(4,4) =  tan(incl*0.5)^fr * cos(omega);
        tm(4,5) =  0.0;
        tm(4,6) =  0.0;

        % ---- partials of psi wrt (a e i O w nu/m/tau)
        tm(5,1) =  0.0;
        tm(5,2) =  0.0;
        tm(5,3) =  cos(omega) * (0.5*tan(incl*0.5)*tan(incl*0.5) + 0.5) * fr * tan(incl*0.5)^(fr-1);
        tm(5,4) =  -tan(incl*0.5)^fr * sin(omega);
        tm(5,5) =  0.0;
        tm(5,6) =  0.0;

        % ---- partials of l wrt (a e i O w nu/m/tau)
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            %[e0,nu] = newtonm ( ecc,anomaly );
            tm(6,1) = 0.0;           
            tm(6,2) = 0.0; 
            tm(6,3) = 0.0;
            tm(6,4) = fr;
            tm(6,5) = 1.0;
            tm(6,6) = 1.0; 
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                tm(6,1) = 0.0;
                tm(6,2) = 0.0;
                tm(6,3) = 0.0;
                tm(6,4) = fr;
                tm(6,5) = 1.0;
                tm(6,6) = 1.0;
            end
        end
        
        % ---------- calculate the output covariance matrix -----------
        [eqcov] =  tm*classcov*tm';

