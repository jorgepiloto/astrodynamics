% ------------------------------------------------------------------------------
%
%                           function EQCM_to_ECI_RTN_sal
%
%  this function finds the interceptor's ECI (RTN) pos/vel vectors 
%  given the ECI target and Modified Equidistant Cylindrical (EQCM) 
%  interceptor vectors.
%
%  Routine IS dependent on km for distance unit
%  all vectors are column vectors
%
%  units are in km, seconds, and radians
%  all vectors are column vectors
%
%  author        : sal alfano         719-573-2600   25 sep 2012
%
% ------------------------------------------------------------------------------

function [r_int_ECI, v_int_ECI] = EQCM_to_ECI_RTN_sal(r_tgt_ECI, v_tgt_ECI, r_int_EQCM, v_int_EQCM )

%  find rotation matrix from ECI to RTN frame
%  convert target and compute position vector magnitude
%  In RTN frame lambda_tgt will be 0
        rot_ECI_to_RTN1 = f_ECI_to_RTN_sal(r_tgt_ECI,v_tgt_ECI);        
        r_tgt_RTN1 = rot_ECI_to_RTN1*r_tgt_ECI;
        v_tgt_RTN1 = rot_ECI_to_RTN1*v_tgt_ECI;
        mag_r_tgt  = norm(r_tgt_RTN1);

%  find necessary orbital elements of target at present (nu1) location 
        mu_m = 398600.4415; % use km!
        h_vec_tgt = cross(r_tgt_RTN1,v_tgt_RTN1);
        p_tgt = dot(h_vec_tgt,h_vec_tgt)/mu_m;
        ecc_vec_tgt = cross(v_tgt_RTN1,h_vec_tgt)/mu_m-r_tgt_RTN1/mag_r_tgt;
        ecc_tgt = norm(ecc_vec_tgt);
        a_tgt= p_tgt/(1-ecc_tgt*ecc_tgt);
        if ecc_tgt > 0.00001
            perigee_unit = ecc_vec_tgt/ecc_tgt;
        else
            perigee_unit = r_tgt_RTN1/mag_r_tgt;
        end;
        lambda_perigee = atan2(perigee_unit(2,1),perigee_unit(1,1));
        nu1 = -lambda_perigee;
        
%  find nu2 and lambda from orbit arc
        arclength = r_int_EQCM(2,1);
        [ea1, m] = newtonnu ( ecc_tgt,nu1 );        
        if abs(arclength) > 0.001   % 1 m tolerance on the arclength
            DE = arclength / a_tgt;  % arclength
            Deltaea = inverselliptic2(DE,ecc_tgt^2 );  % this will be an eccentric anomaly
            [F1, E1] = elliptic12( ea1, ecc_tgt^2 );  
            ea2e = ea1 + Deltaea;
        
            % refine range if ea0 is non-zero because Deltaea is not the same at
            % different points in the orbit
            ii = 1;
            arclength1a = arclength + 10.0;  % initial value
            while (ii < 10) && (abs(arclength1a - arclength) > 0.001)  % 1 m
                [F2, E2] = elliptic12( ea2e, ecc_tgt^2 );   
                arclength1a = a_tgt*(E2 - E1);
                corr = arclength / (ea2e-ea1);   % km/rad
                ea2e = ea2e - (arclength1a-arclength)/corr;
                ii = ii + 1;
            end 
            [m,nu2] = newtone ( ecc_tgt,ea2e );  % convert eccentric back to true   
        else
            nu2 = nu1;   
        end   % if arclength is nearly 0.0
            
        lambda = nu2 - nu1;
        sin_lambda = sin(lambda);
        cos_lambda = cos(lambda);
        
%  find future position and velocity of target using nu2
        r2_tgt = p_tgt / (1.0 + ecc_tgt*cos(nu2));
        P_vec = perigee_unit;
        Q_vec = cross([0 0 1]',P_vec);
        r2_vec_tgt = r2_tgt*(cos(nu2)*P_vec+sin(nu2)*Q_vec);
        v2_vec_tgt = sqrt(mu_m/p_tgt)*(-sin(nu2)*P_vec + (ecc_tgt+cos(nu2))*Q_vec);
        
%  rotate to future target (RTN2) frame
        rot_RTN1_to_RTN2 = f_ECI_to_RTN_sal(r2_vec_tgt,v2_vec_tgt);
        r_tgt_RTN2 = rot_RTN1_to_RTN2*r2_vec_tgt;
        v_tgt_RTN2 = rot_RTN1_to_RTN2*v2_vec_tgt;

%  find phi 
        phi = r_int_EQCM(3,1)/r2_tgt;
        cos_phi = cos(phi);
        sin_phi = sin(phi);
        
%  find interceptor position unit vector in RTN1 frame
        r_int_unit_RTN1(3,1) = sin_phi;
        r_int_unit_RTN1(2,1) = sin_lambda*cos_phi;
        r_int_unit_RTN1(1,1) = cos_lambda*cos_phi;
        
%  find SEZ conversion
        rot_to_SEZ = zeros(3,3);
        rot_to_SEZ(1,1) = sin_phi*cos_lambda;
        rot_to_SEZ(1,2) = sin_phi*sin_lambda;
        rot_to_SEZ(1,3) = -cos_phi;
        rot_to_SEZ(2,1) = -sin_lambda;
        rot_to_SEZ(2,2) = cos_lambda;
        rot_to_SEZ(2,3) = 0.0;
        rot_to_SEZ(3,1) = cos_phi*cos_lambda;
        rot_to_SEZ(3,2) = cos_phi*sin_lambda;
        rot_to_SEZ(3,3) = sin_phi;
        r_int_unit_SEZ = rot_to_SEZ*r_int_unit_RTN1;        

%  determine proper scaling from Z component of interecptor
        r_int_SEZ(3,1) = r_int_EQCM(1,1)+r_tgt_RTN2(1,1);
        mag_r_int = r_int_SEZ(3,1)/r_int_unit_SEZ(3,1);
        r_int_RTN1 = mag_r_int*r_int_unit_RTN1;

%  find velocity component positions in RTN1 frame
        lamda_dot = (v_int_EQCM(2,1) + (v_tgt_RTN1(2,1)/mag_r_tgt)*mag_r_tgt)/r2_tgt;
        v_int_SEZ(1,1) = (-v_int_EQCM(3,1)/r2_tgt)*mag_r_int;
        v_int_SEZ(2,1) = lamda_dot*mag_r_int*cos_phi;
        v_int_SEZ(3,1) = v_int_EQCM(1,1) + v_tgt_RTN2(1,1);
        v_int_RTN1 = rot_to_SEZ'*v_int_SEZ;

%  now rotate all into original frame        
        r_int_ECI = rot_ECI_to_RTN1'*r_int_RTN1;
        v_int_ECI = rot_ECI_to_RTN1'*v_int_RTN1;

