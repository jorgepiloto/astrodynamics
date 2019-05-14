function [ Een, Eex ] = ShadowEntryExit( RSun, rp, a, ecc, incl, raan, argp, nu, mu )

    Een = 0.0;
    Eex = 0.0;
     
    % Semi-Parameter
    p = a*(1.0 - ecc^2);

    % Eccentric Anomaly:
    sinE = (sind(nu)*sqrt(1 - ecc^2))/(1 + ecc*cosd(nu));
    cosE = (ecc + cosd(nu))/(1 + ecc*cosd(nu));

    % (5.3)
    Px = cosd(raan)*cosd(argp) - sind(raan)*sind(argp)*cosd(incl);
    Py = cosd(raan)*sind(argp) + sind(raan)*cosd(argp)*cosd(incl);
    Pz = sind(raan)*sind(incl);
    P_ = [Px,Py,Pz];
    Qx = -sind(raan)*cosd(argp) - cosd(raan)*sind(argp)*cosd(incl);
    Qy = -sind(raan)*sind(argp) + cosd(raan)*cosd(argp)*cosd(incl);
    Qz = cosd(raan)*sind(incl);
    Q_ = [Qx,Qy,Qz];

    % (5.6)
    beta = dot(P_,RSun)/norm(RSun);
    zeta = dot(Q_,RSun)/norm(RSun);
    A0 = ((rp/p)^4)*(ecc^4) - 2*((rp/p)^2)*(zeta^2 - beta^2)*ecc^2 + (beta^2 + zeta^2)^2;
    A1 = 4*((rp/p)^4)*ecc^3 - 4*((rp/p)^2)*(zeta^2 - beta^2)*ecc;
    A2 = 6*((rp/p)^4)*ecc^2 - 2*((rp/p)^2)*(zeta^2 - beta^2) - 2*((rp/p)^2)*(1-zeta^2)*ecc^2 + 2*(zeta^2 - beta^2)*(1 - zeta^2) - 4*(beta^2)*(zeta^2);
    A3 = 4*((rp/p)^4)*ecc - 4*((rp/p)^2)*(1-zeta^2)*ecc;
    A4 = (rp/p)^4 - 2*((rp/p)^2)*(1 - zeta^2) + (1 - zeta^2)^2;

    [r1r,~,r2r,~,r3r,~,r4r,~] = quartic( A0,A1,A2,A3,A4,'R' );
    nu(1) = acosd(r1r);
    nu(2) = acosd(r2r);
    nu(3) = acosd(r3r);
    nu(4) = acosd(r4r);
    k = 1;
    for incl = 1:length(nu)
        check(incl) = beta*cosd(nu(incl)) + zeta*sind(nu(incl));
        if check(incl) < 0
            nugood(k) = nu(incl);
            k = k + 1;
            disp('----------------------------------------------------------------------')
            disp(strcat(['Valid Eccentric Anomaly (beta*cos(nu) + zeta*sin(nu) < 0): ',num2str(nu(incl)),' deg']))
            before = A0*cosd(nugood-0.01)^4 + A1*cosd(nugood-0.01)^3 + A2*cosd(nugood-0.01)^2 + A3*cosd(nugood-0.01) + A4;
            after = A0*cosd(nugood+0.01)^4 + A1*cosd(nugood+0.01)^3 + A2*cosd(nugood+0.01)^2 + A3*cosd(nugood+0.01) + A4;
            if before < 0 && after > 0
                disp('Entering Shadow for This Eccentric Anomaly')
                Een = nu(incl);
            elseif before > 0 && after < 0 
                disp('Exiting Shadow for This Eccentric Anomaly')
                Eex = nu(incl);
            else
                disp('Error - whether entering or exiting is not clear')
            end
        end
    end
end

