% ----------------------------------------------------------------------------
%
%                           function covcl2ct
%
%  this function transforms a six by six covariance matrix expressed in classical elements
%    into one expressed in cartesian elements
%
%  author        : david vallado
%
%  revisions
%    vallado     - simplify code using pqw-eci transformation    12 may 2017
%    vallado     - major update                                  26 aug 2015
%
%  inputs          description                    range / units
%    classcov    - 6x6 classical covariance matrix
%    classstate  - 6x1 classical orbit state      (a e i O w nu/m)
%    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
%
%  outputs       :
%    cartcov     - 6x6 cartesian covariance matrix
%    tm          - transformation matrix
%
%  locals        :
%    r           - matrix of partial derivatives
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omaga       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    m           - mean anomaly                   0.0  to 2pi rad
%    p1,p2,p3,p4 - denominator terms for the partials
%    e0          - eccentric anomaly              0.0  to 2pi rad
%    true1, true2- temp true anomaly              0.0  to 2pi rad
%
%  coupling      :
%    newtonm     - newton iteration for m and ecc to nu
%    newtonnu    - newton iteration for nu and ecc to m
%    constastro
%
%  references    :
%    Vallado and Alfano 2015
%
%   [cartcov,tm] = covcl2ct( classcov,classstate,anom )
% ----------------------------------------------------------------------------

function [cartcov,tm] = covcl2ctnew( classcov, classstate, anom )

       % -------- define gravitational constant
        constastro;

        % --------- determine which set of variables is in use ---------
        % ---- parse the input vector into the classical elements -----
        a = classstate(1);
        n = sqrt(mum/a^3);
        ecc     = classstate(2);
        incl    = classstate(3);
        omega   = classstate(4);
        argp    = classstate(5);
        % -------- if mean anomaly is used, convert to true anomaly
        % -------- eccentric anomaly (e) is needed for both
        if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
            mean = classstate(6);
           [e, nu] = newtonm(ecc, mean);
        else
            if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            % note that mean is not used in the partials, but nu is! 
               nu = classstate(6);
               [e, mean] = newtonnu(ecc, nu);
            end
        end
        
        p = a*(1-ecc^2)/1000;  % needs to be in km
        [r,v] = coe2rv ( p,ecc,incl,omega,argp,nu,0.0,0.0,0.0 );
        rx = r(1)*1000;
        ry = r(2)*1000;
        rz = r(3)*1000;
        vx = v(1)*1000;
        vy = v(2)*1000;
        vz = v(3)*1000;
        
        % assign trig values for efficiency
        sin_inc = sin(incl);
        cos_inc = cos(incl);
        sin_anode = sin(omega);
        cos_anode = cos(omega);
        sin_w = sin(argp);
        cos_w = cos(argp);
        sin_nu = sin(nu);
        cos_nu = cos(nu);
       
        % assign elements of PQW to ECI transformation (pg 168)
        p11 =  cos_anode*cos_w - sin_anode*sin_w*cos_inc;
        p12 = -cos_anode*sin_w - sin_anode*cos_w*cos_inc;
        p13 =  sin_anode*sin_inc;
        p21 =  sin_anode*cos_w + cos_anode*sin_w*cos_inc;
        p22 = -sin_anode*sin_w + cos_anode*cos_w*cos_inc;
        p23 = -cos_anode*sin_inc;
        p31 =  sin_w*sin_inc;
        p32 =  cos_w*sin_inc;
        p33 =  cos_inc;

        % assign constants for efficiency
        p0 =  sqrt(mum / (a*(1.0 - ecc*ecc)));
        p1 = (1.0 - ecc*ecc) / (1.0 + ecc*cos_nu);
        p2 = 1.0 / (2.0*a) * p0;
        p3 = (2.0*a*ecc + a*cos_nu + a*cos_nu*ecc*ecc) / ((1.0 + ecc*cos_nu)^2);
        p4 = ecc*mum / (a*(1.0 - ecc*ecc)^2 *p0);
        p5 = a * p1;
        p6 = a*(1.0 - ecc*ecc) / ((1.0 + ecc*cos_nu)^2);
        
        p8 = (1.0 - ecc^2)^1.5 / ( (1.0 + ecc*cos(nu))^2 );  % dm/dv
        p9 = -sin(nu)*((ecc*cos(nu) + 1)*(ecc+cos(nu))/sqrt((ecc + cos(nu))^2) + 1.0 - 2.0*ecc^2 - ecc^3*cos(nu)) / ( (ecc*cos(nu) + 1.0)^2 * sqrt(1-ecc^2) );  % dm/de   
        

        % ---------------- calculate matrix elements ------------------
        % ---- partials of rx wrt (a e i O w m)
        tm(1,1) = p1 * (p11*cos_nu + p12*sin_nu);
        tm(1,1) = rx/a; % alternate approach
        tm(1,2) = -p3 * (p11*cos_nu + p12*sin_nu);
        tm(1,3) = p5 * p13*(sin_w*cos_nu + cos_w*sin_nu);
        tm(1,4) = -p5 * (p21*cos_nu + p22*sin_nu);
        tm(1,5) = p5 * (p12*cos_nu - p11*sin_nu);
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            tm(1,6) = p6 * (-p11*sin_nu + p12*(ecc + cos_nu))
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p10 = a * (ecc^2 - 1.0) / (ecc*cos(nu) + 1.0)^2;
                tm(1,6) = p10 * (ecc*cos(omega)*sin(argp)+cos(omega)*cos(argp)*sin(nu)+cos(omega)*sin(argp)*cos(nu)+ecc*cos(incl)*sin(omega)*cos(argp)+cos(incl)*sin(omega)*cos(argp)*cos(nu)-cos(incl)*sin(omega)*sin(argp)*sin(nu));
                tm(1,6) = tm(1,6) / p8;
                tm(1,2) = tm(1,2) - tm(1,6) * p9;
            end
        end
        
        % ---- partials of ry wrt (a e i O w nu/m)
        tm(2,1) = p1 * (p21*cos_nu + p22*sin_nu);
        tm(2,1) = ry/a; 
        tm(2,2) = -p3 * (p21*cos_nu + p22*sin_nu);
        tm(2,3) = p5 * p23*(sin_w*cos_nu + cos_w*sin_nu);
        tm(2,4) = p5 * (p11*cos_nu + p12*sin_nu);
        tm(2,5) = p5 * (p22*cos_nu - p21*sin_nu);
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            tm(2,6) = p6 * (-p21*sin_nu + p22*(ecc + cos_nu))
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p10 = a * (ecc^2 - 1.0) / (ecc*cos(nu) + 1.0)^2;
                tm(2,6) = p10 * (ecc*sin(omega)*sin(argp)+ sin(omega)*cos(argp)*sin(nu)+ sin(omega)*sin(argp)*cos(nu)- ecc*cos(incl)*cos(omega)*cos(argp)- cos(incl)*cos(omega)*cos(argp)*cos(nu)+ cos(incl)*cos(omega)*sin(argp)*sin(nu));
                tm(2,6) = tm(2,6) / p8;
                tm(2,2) = tm(2,2) - tm(2,6) * p9;
           end
        end
               
        % ---- partials of rz wrt (a e i O w nu/m)
        tm(3,1) = p1 * (p31*cos_nu + p32*sin_nu);
        tm(3,1) = rz/a; 
        tm(3,2) = -p3 * sin_inc * (cos_w*sin_nu + sin_w*cos_nu);
        tm(3,3) =  p5 * cos_inc * (cos_w*sin_nu + sin_w*cos_nu);
        tm(3,4) = 0.0;
        tm(3,5) = p5 * sin_inc*(cos_w*cos_nu - sin_w*sin_nu);
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            tm(3,6) = p6 * (-p31*sin_nu + p32*(ecc + cos_nu))
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p10 = -a * (ecc^2 - 1.0) / (ecc*cos(nu) + 1.0)^2;
                tm(3,6) = p10 * sin(incl)*(cos(argp+nu)+ecc*cos(argp));
                tm(3,6) = tm(3,6) / p8;
                tm(3,2) = tm(3,2) - tm(3,6) * p9;
            end
        end
                
        % ---- partials of vx wrt (a e i O w nu/m)
        tm(4,1) = p2 * (p11*sin_nu - p12*(ecc + cos_nu));
        tm(4,1) = -vx/(2.0*a);  
        tm(4,2) = -p4 * (p11*sin_nu - p12*(ecc + cos_nu)) + p12*p0;
        tm(4,3) = -p0 * sin_anode*(p31*sin_nu - p32*(ecc + cos_nu));
        tm(4,4) = p0 * (p21*sin_nu - p22*(ecc + cos_nu));
        tm(4,5) = -p0 * (p12*sin_nu + p11*(ecc + cos_nu));
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            tm(4,6) = -p0 * (p11*cos_nu + p12*sin_nu)
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p10 = sqrt(-mum/(a*(ecc^2-1.0)));
                tm(4,6) = p10 * (cos(omega)*sin(argp)*sin(nu)-cos(omega)*cos(argp)*cos(nu)+cos(incl)*sin(omega)*cos(argp)*sin(nu)+cos(incl)*sin(omega)*sin(argp)*cos(nu));
                tm(4,6) = tm(4,6) / p8;
                tm(4,2) = tm(4,2) - tm(4,6) * p9;
            end
        end
        
        % ---- partials of vy wrt (a e i O w nu/m)
        tm(5,1) = p2 * (p21*sin_nu - p22*(ecc + cos_nu));
        tm(5,1) = -vy/(2.0*a);  
        tm(5,2) = -p4 * (p21*sin_nu - p22*(ecc + cos_nu)) + p22*p0;
        tm(5,3) = p0 * cos_anode*(p31*sin_nu - p32*(ecc + cos_nu));
        tm(5,4) = p0 * (-p11*sin_nu + p12*(ecc + cos_nu));
        tm(5,5) = -p0 * (p22*sin_nu + p21*(ecc + cos_nu));
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            tm(5,6) = -p0 * (p21*cos_nu + p22*sin_nu)
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p10 = sqrt(-mum/(a*(ecc^2-1.0)));
                tm(5,6) = -p10 * (sin(omega)*cos(argp)*cos(nu)-sin(omega)*sin(argp)*sin(nu)+cos(incl)*cos(omega)*cos(argp)*sin(nu)+cos(incl)*cos(omega)*sin(argp)*cos(nu));
                tm(5,6) = tm(5,6) / p8;
                tm(5,2) = tm(5,2) - tm(5,6) * p9;
            end
        end
        
        % ---- partials of vz wrt (a e i O w nu/m)
        tm(6,1) = p2 * (p31*sin_nu - p32*(ecc + cos_nu));
        tm(6,1) = -vz/(2.0*a);  
        tm(6,2) = -p4 * (p31*sin_nu - p32*(ecc + cos_nu)) + p32*p0;
        tm(6,3) = p0 * cos_inc*(cos_w*cos_nu - sin_w*sin_nu + ecc*cos_w);
        tm(6,4) = 0.0;
        tm(6,5) = -p0 * (p32*sin_nu + p31*(ecc + cos_nu));
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            tm(6,6) = -p0 * (p31*cos_nu + p32*sin_nu)
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p10 = sqrt(-mum/(a*(ecc^2-1.0)));
                tm(6,6) = p10 * (-sin(incl)*sin(argp+nu));
                tm(6,6) = tm(6,6) / p8;
                tm(6,2) = tm(6,2) - tm(6,6) * p9;
            end
        end
        
        % ---------- calculate the output covariance matrix -----------
        cartcov =  tm*classcov*tm';

