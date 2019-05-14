% ----------------------------------------------------------------------------
%
%                           function covct2eq
%
%  this function transforms a six by six covariance matrix expressed in
%    cartesian vectors into one expressed in equinoctial elements.
%
%  author        : david vallado                  719-573-2600   14 jul 2002
%
%  revisions
%    vallado     - major update                                  26 aug 2015
%
%  inputs          description                    range / units
%    cartcov     - 6x6 cartesian covariance matrix
%    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
%    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
%    fr          - retrograde factor               +1, -1
%
%  outputs       :
%    eqcov       - 6x6 equinoctial covariance matrix
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
%    e0          - eccentric anomaly              0.0  to 2pi rad
%    tau         - time from perigee passage
%    n           - mean motion                    rad
%    af          - component of ecc vector
%    ag          - component of ecc vector
%    chi         - component of node vector in eqw
%    psi         - component of node vector in eqw
%    meanlon     - mean longitude                 rad
%
%  coupling      :
%    constastro
%
%  references    :
%    Vallado and Alfano 2015
%
% [eqcov, tm] = covct2eq ( cartcov, cartstate, anom, fr);
% ----------------------------------------------------------------------------
 
function [eqcov, tm] = covct2eq ( cartcov, cartstate, anom, fr)

        % -------- define gravitational constant
        constastro;

        % -------- parse the input vectors into cartesian and classical components
        rx = cartstate(1) * 1000.0;
        ry = cartstate(2) * 1000.0;
        rz = cartstate(3) * 1000.0;
        vx = cartstate(4) * 1000.0;
        vy = cartstate(5) * 1000.0;
        vz = cartstate(6) * 1000.0;
        reci = [rx; ry; rz];
        veci = [vx; vy; vz];
        magr = mag(reci); % m
        magv = mag(veci); % m
               
        a = 1.0 / (2.0/magr - magv^2/mum);
        n = sqrt(mum / a^3);

        hx =  ry*vz - rz*vy;
        hy = -rx*vz + rz*vx;
        hz =  rx*vy - ry*vx;
        h_vec = [hx hy hz]';
        %h_vec = cross(reci, veci)  %same
        w_vec = h_vec/mag(h_vec);
        chi = w_vec(1)/(1.0 + fr * w_vec(3));
        psi = -w_vec(2)/(1.0 + fr * w_vec(3));

        % components of equinoctial system
        p0 = 1.0/(1.0 + chi^2 + psi^2);
        fe = p0 * (1.0 - chi^2 + psi^2);  
        fq = p0 * 2.0 * chi * psi;
        fw = p0 * -2.0 * fr * chi;
        ge = p0 * 2.0 * fr* chi * psi;
        gq = p0 * fr * (1.0 + chi^2 - psi^2);
        gw = p0 * 2.0 * psi;
        f_vec = [fe fq fw]';
        g_vec = [ge gq gw]';
        we = w_vec(1);
        wq = w_vec(2);
        ww = w_vec(3);

        r_dot_v = dot(reci, veci); 
        p1 = magv*magv - mum/magr;
        ecc_x = (p1*rx - r_dot_v*vx)/mum;
        ecc_y = (p1*ry - r_dot_v*vy)/mum;
        ecc_z = (p1*rz - r_dot_v*vz)/mum;
        ecc_vec = [ecc_x ecc_y ecc_z]';

        af = dot(ecc_vec, f_vec);
        ag = dot(ecc_vec, g_vec);
        
        X = dot(reci, f_vec);
        Y = dot(reci, g_vec);

        b = 1.0 / (1.0 + sqrt(1.0 - af^2 - ag^2));
        p0 = 1.0 / (a*sqrt(1.0 - af^2 - ag^2));
        sinF = ag + p0*((1.0 - ag^2*b)*Y - ag*af*b*X);  
        cosF = af + p0*((1.0 - af^2*b)*X - ag*af*b*Y);          
        F = atan2(sinF,cosF);
        if F < 0.0 
            F = F + 2.0 * pi;
        end

        XD = n*a^2/magr * ( af*ag*b*cos(F) - (1.0 - ag^2*b)*sin(F) );
        YD = n*a^2/magr * ( (1.0 - af^2*b)*cos(F) - af*ag*b*sin(F) );
        
        A = sqrt(mum*a);
        B = sqrt(1.0 - ag^2 - af^2);
        C = 1.0 + chi^2 + psi^2;
       
        partXDaf =  a*XD*YD / (A*B) - A/(magr^3)*( a*ag*X/(1 + B) + X*Y/B ); 
        partYDaf = -a*XD^2 / (A*B) - A/(magr^3)*( a*ag*Y/(1 + B) - X^2/B ); 
        partXDag =  a*YD^2 / (A*B) + A/(magr^3)*( a*af*X/(1 + B) - Y^2/B );
        partYDag = -a*XD*YD / (A*B) + A/(magr^3)*( a*af*Y/(1 + B) + X*Y/B ); 
         
        % ---------------- calculate matrix elements ------------------
        % ---- partials of a wrt (rx ry rz vx vy vz)
        if strcmp(anom,'truea') == 1 || strcmp(anom,'meana') == 1  % 1 is true
            p0 = 2.0*a^2 / magr^3;
            p1 = 2.0 / (n^2*a);
        else
            if strcmp(anom,'truen') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p0 = -3.0*n*a/magr^3; 
                p1 = -3.0 / (n*a^2); 
            end
        end
        
        tm(1,1) = p0*rx;
        tm(1,2) = p0*ry;
        tm(1,3) = p0*rz;
        tm(1,4) = p1*vx; 
        tm(1,5) = p1*vy;
        tm(1,6) = p1*vz;

        % ---- partials of v wrt ag
        tm34 = partXDag*fe + partYDag*ge;
        tm35 = partXDag*fq + partYDag*gq;
        tm36 = partXDag*fw + partYDag*gw;      
        % ---- partials of af wrt (rx ry rz vx vy vz)
        p0 = 1.0 / mum;
        tm(2,1) = -a*b*af*B*rx/(magr^3) - (ag*(chi*XD - psi*fr*YD)*we)/(A*B) + (B/A)*tm34;
        tm(2,2) = -a*b*af*B*ry/(magr^3) - (ag*(chi*XD - psi*fr*YD)*wq)/(A*B) + (B/A)*tm35;
        tm(2,3) = -a*b*af*B*rz/(magr^3) - (ag*(chi*XD - psi*fr*YD)*ww)/(A*B) + (B/A)*tm36;
        tm(2,4) = p0*((2.0*X*YD - XD*Y)*ge - Y*YD*fe)  - (ag*(psi*fr*Y - chi*X)*we) / (A*B);
        tm(2,5) = p0*((2.0*X*YD - XD*Y)*gq - Y*YD*fq)  - (ag*(psi*fr*Y - chi*X)*wq) / (A*B);
        tm(2,6) = p0*((2.0*X*YD - XD*Y)*gw - Y*YD*fw)  - (ag*(psi*fr*Y - chi*X)*ww) / (A*B);

        % ---- partials of v wrt af
        tm24 = partXDaf*fe + partYDaf*ge;
        tm25 = partXDaf*fq + partYDaf*gq;
        tm26 = partXDaf*fw + partYDaf*gw;
        % ---- partials of ag wrt (rx ry rz vx vy vz)
        p0 = 1.0 / mum;
        tm(3,1) = -a*b*ag*B*rx/(magr^3) + (af*(chi*XD - psi*fr*YD)*we)/(A*B) - (B/A)*tm24;
        tm(3,2) = -a*b*ag*B*ry/(magr^3) + (af*(chi*XD - psi*fr*YD)*wq)/(A*B) - (B/A)*tm25;
        tm(3,3) = -a*b*ag*B*rz/(magr^3) + (af*(chi*XD - psi*fr*YD)*ww)/(A*B) - (B/A)*tm26;
        tm(3,4) = p0*((2.0*XD*Y - X*YD)*fe - X*XD*ge)  + (af*(psi*fr*Y - chi*X)*we) / (A*B);
        tm(3,5) = p0*((2.0*XD*Y - X*YD)*fq - X*XD*gq)  + (af*(psi*fr*Y - chi*X)*wq) / (A*B);
        tm(3,6) = p0*((2.0*XD*Y - X*YD)*fw - X*XD*gw)  + (af*(psi*fr*Y - chi*X)*ww) / (A*B);
       
        % ---- partials of chi wrt (rx ry rz vx vy vz)
        tm(4,1) = -C*YD*we / (2.0*A*B);
        tm(4,2) = -C*YD*wq / (2.0*A*B);
        tm(4,3) = -C*YD*ww / (2.0*A*B);
        tm(4,4) = C*Y*we / (2.0*A*B);
        tm(4,5) = C*Y*wq / (2.0*A*B);
        tm(4,6) = C*Y*ww / (2.0*A*B);

        % ---- partials of psi wrt (rx ry rz vx vy vz)
        tm(5,1) = -fr*C*XD*we / (2.0*A*B);
        tm(5,2) = -fr*C*XD*wq / (2.0*A*B);
        tm(5,3) = -fr*C*XD*ww / (2.0*A*B);
        tm(5,4) = fr*C*X*we / (2.0*A*B);
        tm(5,5) = fr*C*X*wq / (2.0*A*B);
        tm(5,6) = fr*C*X*ww / (2.0*A*B);

        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
% not ready yet
%            p0 = -sign(argp)/sqrt(1-cos(argp)*cos(argp));
%            p1 = 1.0 / mum;
%             tm(6,1) = p0*(p1*()/(n*ecc) - ()/n*()/(n^2*ecc) - tm(ecc/ry*()) + fr*-vz*nodey/n^2 + ...
%                      ;
%             tm(6,2) = p0*();
%             tm(6,3) = p0*();
%             tm(6,4) = p0*();
%             tm(6,5) = p0*();
%             tm(6,6) = p0*();
             tm(6,1) = 0;
             tm(6,2) = 0;
             tm(6,3) = 0;
             tm(6,4) = 0;
             tm(6,5) = 0;
             tm(6,6) = 0;
         else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                % ---- partials of meanlon wrt (rx ry rz vx vy vz)
                tm(6,1) = -vx/A + (chi*XD - psi*fr*YD)*we/(A*B) - (b*B/A)*(ag*tm34 + af*tm24);
                tm(6,2) = -vy/A + (chi*XD - psi*fr*YD)*wq/(A*B) - (b*B/A)*(ag*tm35 + af*tm25);
                tm(6,3) = -vz/A + (chi*XD - psi*fr*YD)*ww/(A*B) - (b*B/A)*(ag*tm36 + af*tm26);
                tm(6,4) = -2.0*rx/A + (af*tm(3,4) - ag*tm(2,4))/(1.0 + B) + (fr*psi*Y - chi*X)*we/A;
                tm(6,5) = -2.0*ry/A + (af*tm(3,5) - ag*tm(2,5))/(1.0 + B) + (fr*psi*Y - chi*X)*wq/A;
                tm(6,6) = -2.0*rz/A + (af*tm(3,6) - ag*tm(2,6))/(1.0 + B) + (fr*psi*Y - chi*X)*ww/A;
            end
        end

        % ---------- calculate the output covariance matrix -----------
        [eqcov] =  tm*cartcov*tm';
        


