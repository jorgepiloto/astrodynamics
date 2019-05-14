% ----------------------------------------------------------------------------
%
%                           function covct2cl
%
%  this function transforms a six by six covariance matrix expressed in cartesian elements
%    into one expressed in classical elements
%
%  author        : david vallado                  719-573-2600   21 jun 2002
%
%  revisions
%    vallado     - major update                                  26 aug 2015
%
%  inputs          description                    range / units
%    cartcov     - 6x6 cartesian covariance matrix
%    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
%    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
%
%  outputs       :
%    classcov    - 6x6 classical covariance matrix
%    tm          - transformation matrix
%
%  locals        :
%    r           - matrix of partial derivatives
%    rj2000      - position vector                km
%    x,y,z       - components of position vector  km
%    vj2000      - velocity vector                km/s
%    vx,vy,vz    - components of position vector  km/s
%    p           - semilatus rectum               km
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omaga       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    m           - mean anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
%    magr        - magnitude of position vector   km
%    magv        - magnitude of velocity vector   km/s
%
%  coupling      :
%    constastro
%    rv2coe      - position and velocity vectors to classical elements
%
%  references    :
%    Vallado and Alfano 2015
%
%   [classcov,tm] = covct2cl( cartcov,cartstate,anom );
% ----------------------------------------------------------------------------

function [classcov,tm] = covct2clnew(cartcov, cartstate, anom )

        % -------- define gravitational constant
        constastro;

        % -------- parse the input vectors into cartesian and classical components
        rx  = cartstate(1) * 1000.0;
        ry  = cartstate(2) * 1000.0;
        rz  = cartstate(3) * 1000.0;
        vx = cartstate(4) * 1000.0;
        vy = cartstate(5) * 1000.0;
        vz = cartstate(6) * 1000.0;
        reci = [rx/1000;ry/1000;rz/1000];
        veci = [vx/1000;vy/1000;vz/1000];

        % -------- convert to a classical orbit state for ease of computation
        [p, a, ecc, incl ,omega, argp, nu, mean, arglat, truelon, lonper ] = rv2coe  (reci, veci);
        p = p * 1000.0;
        a = a * 1000.0;
        n = sqrt(mum/a^3);

        % -------- calculate common quantities
        sqrt1me2 = sqrt(1.0 - ecc*ecc);
        magr = sqrt(rx^2 + ry^2 + rz^2);
        magr3 = magr^3;
        magv = sqrt(vx^2 + vy^2 + vz^2);

        % ----------  form pqw position and velocity vectors ----------  
        r_dot_v = dot(reci, veci)*1000*1000;
        
        ecc_term = magv*magv - mum/magr;
        ecc_x = (ecc_term*rx - r_dot_v*vx)/mum;
        ecc_y = (ecc_term*ry - r_dot_v*vy)/mum;
        ecc_z = (ecc_term*rz - r_dot_v*vz)/mum;
        ecc_vec = [ecc_x ecc_y ecc_z]';
        
        hx = ry*vz - rz*vy;
        hy = rz*vx - rx*vz;
        hz = rx*vy - ry*vx;
        h_vec = [hx hy hz]';
        h = mag(h_vec);
        h_squared = h*h;
        
        nx = -hy;
        ny = hx;
        nz = 0.0;
        node_vec = [nx ny nz]';
        node = mag(node_vec);
        n_squared = node*node;
        n_dot_e = dot(node_vec,ecc_vec);
              
        sign_anode = sign(ny);
        cos_anode = nx/node;
        omega = sign_anode*acos(cos_anode);
        
        sign_w = sign((magv^2 - mum/magr)*rz-r_dot_v*vz);
        cos_w = n_dot_e/(ecc*node);
        argp = sign_w*acos(cos_w);
        w_scale = -sign_w/sqrt(1-cos_w*cos_w);
        
        r_dot_e = dot(reci,ecc_vec)*1000;
        cos_nu = r_dot_e/(magr*ecc);
        sign_nu = sign(r_dot_v);
        nu = sign_nu*acos(cos_nu);
        nu_scale = -sign_nu/sqrt(1-cos_nu*cos_nu);
        
        % ---------------- calculate matrix elements ------------------
        % ---- partials of a wrt (x y z vx vy vz)
        p0 = 2.0*a^2 / magr^3;
        p1 = 2.0 / (n^2*a);
        tm(1,1) = p0*rx;
        tm(1,2) = p0*ry;
        tm(1,3) = p0*rz;
        tm(1,4) = p1*vx; 
        tm(1,5) = p1*vy;
        tm(1,6) = p1*vz;
       
        % ---- partials of ecc wrt (x y z vx vy vz)
        p0 = 1.0 / (mum*ecc);
        tm(2,1) = -p0*(((vx*vy - mum*rx*ry/magr3)*ecc_y) + ((vx*vz - mum*rx*rz/magr3)*ecc_z) - (vy*vy + vz*vz - mum/magr + mum*rx*rx/magr3)*ecc_x);
        tm(2,2) = -p0*(((vx*vy - mum*rx*ry/magr3)*ecc_x) + ((vy*vz - mum*ry*rz/magr3)*ecc_z) - (vx*vx + vz*vz - mum/magr + mum*ry*ry/magr3)*ecc_y);
        tm(2,3) = -p0*(((vx*vz - mum*rx*rz/magr3)*ecc_x) + ((vy*vz - mum*ry*rz/magr3)*ecc_y) - (vy*vy + vx*vx - mum/magr + mum*rz*rz/magr3)*ecc_z);
        tm(2,4) = -p0*((rx*vy - 2*ry*vx)*ecc_y + (ry*vy + rz*vz)*ecc_x + (rx*vz - 2*rz*vx)*ecc_z);
        tm(2,5) = -p0*((ry*vx - 2*rx*vy)*ecc_x + (rx*vx + rz*vz)*ecc_y + (ry*vz - 2*rz*vy)*ecc_z);
        tm(2,6) = -p0*((rx*vx + ry*vy)*ecc_z   + (rz*vx - 2*rx*vz)*ecc_x + (rz*vy - 2*ry*vz)*ecc_y);
         
        % ---- partials of incl wrt (x y z vx vy vz)
        p3 = 1.0 / node;
        tm(3,1) = -p3*(vy - hz*(vy*hz - vz*hy)/h_squared);
        tm(3,2) =  p3*(vx - hz*(vx*hz - vz*hx)/h_squared);
        tm(3,3) = -p3*(hz*(vy*hx - vx*hy)/h_squared);
        tm(3,4) =  p3*(ry - hz*(ry*hz - rz*hy)/h_squared);
        tm(3,5) = -p3*(rx - hz*(rx*hz - rz*hx)/h_squared);
        tm(3,6) =  p3*(hz*(ry*hx - rx*hy)/h_squared);
   
        % ---- partials of node wrt (x y z vx vy vz)
        p4 = 1.0 / n_squared; 
        tm(4,1) = -p4*vz*ny;
        tm(4,2) =  p4*vz*nx;
        tm(4,3) =  p4*(vx*ny - vy*nx);
        tm(4,4) =  p4*rz*ny;
        tm(4,5) = -p4*rz*nx;
        tm(4,6) =  p4*(ry*nx - rx*ny);
        
        % ---- partials of argp wrt (x y z vx vy vz)
        p5 = 1.0 / (node * a * a);
        temp = -hy*(vy*vy + vz*vz - mum/magr + mum*rx*rx/magr3);
        temp = temp-hx*(vx*vy - mum*rx*ry/magr3) + vz*mum*ecc_x;
        temp = temp/(mum*node*ecc) + vz*hy*n_dot_e/(node*node*node*ecc) - tm(2,1)*n_dot_e/(node*ecc*ecc);
        tm(5,1) = temp*w_scale;
        temp = hx*(vx*vx + vz*vz - mum/magr + mum*ry*ry/magr3);
        temp = temp+hy*(vx*vy - mum*rx*ry/magr3) + vz*mum*ecc_y;
        temp = temp/(mum*node*ecc) - vz*hx*n_dot_e/(node*node*node*ecc) - tm(2,2)*n_dot_e/(node*ecc*ecc);
        tm(5,2) = temp*w_scale;
        temp = -hy*(vx*vz - mum*rx*rz/magr3) + hx*(vy*vz - mum*ry*rz/magr3) + vx*mum*ecc_x + vy*mum*ecc_y;
        temp = -temp/(mum*node*ecc) + (vy*hx - vx*hy)*n_dot_e/(node*node*node*ecc) - tm(2,3)*n_dot_e/(node*ecc*ecc);
        tm(5,3) = temp*w_scale;
        temp = (rx*vy - 2*ry*vx)*hx - hy*(ry*vy + rz*vz) + rz*mum*ecc_x;
        temp = -temp/(mum*node*ecc) - rz*hy*n_dot_e/(node*node*node*ecc) - tm(2,4)*n_dot_e/(node*ecc*ecc);
        tm(5,4) = temp*w_scale;
        temp = -(ry*vx - 2*rx*vy)*hy + hx*(rx*vx + rz*vz) + rz*mum*ecc_y;
        temp = -temp/(mum*node*ecc) + rz*hx*n_dot_e/(node*node*node*ecc) - tm(2,5)*n_dot_e/(node*ecc*ecc);
        tm(5,5) = temp*w_scale;
        temp = -(rz*vx - 2*rx*vz)*hy + hx*(rz*vy - 2*ry*vz) - rx*mum*ecc_x - ry*mum*ecc_y;
        temp = -temp/(mum*node*ecc) + (rx*hy - ry*hx)*n_dot_e/(node*node*node*ecc) - tm(2,6)*n_dot_e/(node*ecc*ecc);
        tm(5,6) = temp*w_scale;    
        
        % ---- partials of true anomaly wrt (x y z vx vy vz)
        temp = ry*(vx*vy - mum*rx*ry/magr3) - rx*ecc_term + rz*(vx*vz - mum*rx*rz/magr3);
        temp = temp - rx*(vy*vy + vz*vz - mum/magr + mum*rx*rx/magr3) + vx*r_dot_v;
        temp = -temp/(mum*magr*ecc) - rx*r_dot_e/(magr3*ecc) - tm(2,1)*r_dot_e/(magr*ecc*ecc);
        tm(6,1) = temp*nu_scale;
        temp = rx*(vx*vy - mum*rx*ry/magr3) - ry*ecc_term + rz*(vy*vz - mum*ry*rz/magr3);
        temp = temp-ry*(vx*vx + vz*vz - mum/magr + mum*ry*ry/magr3) + vy*r_dot_v;
        temp = -temp/(mum*magr*ecc) - ry*r_dot_e/(magr3*ecc) - tm(2,2)*r_dot_e/(magr*ecc*ecc);
        tm(6,2) = temp*nu_scale;
        temp = rx*(vx*vz - mum*rx*rz/magr3) - rz*ecc_term + ry*(vy*vz - mum*ry*rz/magr3);
        temp = temp - rz*(vx*vx + vy*vy - mum/magr + mum*rz*rz/magr3) + vz*r_dot_v;
        temp = -temp/(mum*magr*ecc) - rz*r_dot_e/(magr3*ecc) - tm(2,3)*r_dot_e/(magr*ecc*ecc);
        tm(6,3) = temp*nu_scale;
        temp = ry*(rx*vy - 2*ry*vx) + rx*(ry*vy + rz*vz) + rz*(rx*vz - 2*rz*vx);
        temp = -temp/(mum*magr*ecc) - tm(2,4)*r_dot_e/(magr*ecc*ecc);
        tm(6,4) = temp*nu_scale;
        temp = rx*(ry*vx - 2*rx*vy) + ry*(rx*vx + rz*vz) + rz*(ry*vz - 2*rz*vy);
        temp = -temp/(mum*magr*ecc) - tm(2,5)*r_dot_e/(magr*ecc*ecc);
        tm(6,5) = temp*nu_scale;
        temp = rz*(rx*vx + ry*vy) + rx*(rz*vx - 2*rx*vz) + ry*(rz*vy - 2*ry*vz);
        temp = -temp/(mum*magr*ecc) - tm(2,6)*r_dot_e/(magr*ecc*ecc);
        tm(6,6) = temp*nu_scale;

%       % same answers as above
%       % ---- partials of true anomaly wrt (x y z vx vy vz)
%       p8 = magr^2*magv^2 - mum*magr - r_dot_v^2;
%       p9 = 1.0/(p8^2 + r_dot_v^2 * h^2);
%       tm(6,1) = p9 * ( p8 * (h*vx + r_dot_v*(vy*hz - vz*hy)/h) - r_dot_v*h*(2*rx*magv^2 - mum*rx/magr - 2*r_dot_v*vx) );
%       tm(6,2) = p9 * ( p8 * (h*vy + r_dot_v*(vz*hx - vx*hz)/h) - r_dot_v*h*(2*ry*magv^2 - mum*ry/magr - 2*r_dot_v*vy) );
%       tm(6,3) = p9 * ( p8 * (h*vz + r_dot_v*(vx*hy - vy*hx)/h) - r_dot_v*h*(2*rz*magv^2 - mum*rz/magr - 2*r_dot_v*vz) );
%       tm(6,4) = p9 * ( p8 * (h*rx + r_dot_v*(rz*hy - ry*hz)/h) - r_dot_v*h*(2*vx*magr^2 - 2*r_dot_v*rx) );
%       tm(6,5) = p9 * ( p8 * (h*ry + r_dot_v*(rx*hz - rz*hx)/h) - r_dot_v*h*(2*vy*magr^2 - 2*r_dot_v*ry) );
%       tm(6,6) = p9 * ( p8 * (h*rz + r_dot_v*(ry*hx - rx*hy)/h) - r_dot_v*h*(2*vz*magr^2 - 2*r_dot_v*rz) );

        if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
            % ---- partials of mean anomaly wrt (x y z vx vy vz)
            % then update for mean anomaly
            ecc = mag(ecc_vec);
            p5 = (1.0 - ecc^2)^1.5 / ( (1.0 + ecc*cos(nu))^2 );  % dm/dv
            p6 = -sin(nu)*((ecc*cos(nu) + 1)*(ecc+cos(nu))/sqrt((ecc + cos(nu))^2) + 1.0 - 2.0*ecc^2 - ecc^3*cos(nu)) / ( (ecc*cos(nu) + 1.0)^2 * sqrt(1-ecc^2) );  % dm/de
            tm(6,1) = tm(6,1)*p5 + tm(2,1)*p6;
            tm(6,2) = tm(6,2)*p5 + tm(2,2)*p6;
            tm(6,3) = tm(6,3)*p5 + tm(2,3)*p6;
            tm(6,4) = tm(6,4)*p5 + tm(2,4)*p6;
            tm(6,5) = tm(6,5)*p5 + tm(2,5)*p6;
            tm(6,6) = tm(6,6)*p5 + tm(2,6)*p6;
        end 
        
        % ---------- calculate the output covariance matrix -----------
        classcov = tm*cartcov*tm';

        

