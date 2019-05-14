% ------------------------------------------------------------------------------
%
%                           function ECI_to_EQCM_RTN_sal
%
%  this function finds the relative pos/vel vectors in the Modified 
%  Equidistant Cylindrical (EQCM) frame given the ECI target and 
%  interceptor vectors with RTN (=RSW) ordered components.
%
%  Routine IS dependent on km for distance unit due to mu
%  all vectors are column vectors
%
%  units are in meters, seconds, and radians
%  all vectors are column vectors
%
%  author        : sal alfano              719-573-2600   25 sep 2012
% ------------------------------------------------------------------------------

function [r_int_EQCM, v_int_EQCM] = ECI_to_EQCM_RTN_sal(r_tgt_ECI, v_tgt_ECI, r_int_ECI, v_int_ECI, outfilehill)

%  find rotation matrix from ECI to RTN1 frame for target
%  convert target and interceptor, compute vector magnitudes
        rot_ECI_to_RTN1 = f_ECI_to_RTN_sal(r_tgt_ECI,v_tgt_ECI);        
        r_tgt_RTN1 = rot_ECI_to_RTN1*r_tgt_ECI;
        v_tgt_RTN1 = rot_ECI_to_RTN1*v_tgt_ECI;
        r_int_RTN1 = rot_ECI_to_RTN1*r_int_ECI;
        v_int_RTN1 = rot_ECI_to_RTN1*v_int_ECI;
        mag_r_tgt = norm(r_tgt_RTN1);
        mag_r_int = norm(r_int_RTN1);

%  find lamda/phi (long/lat) rotation angles
%  to go from target to interceptor (lambda_tgt will be 0)
        sin_phi = r_int_RTN1(3,1) / mag_r_int;
        phi = asin(sin_phi);
        cos_phi = cos(phi);
        lambda = atan2(r_int_RTN1(2,1),r_int_RTN1(1,1));
        cos_lambda = cos(lambda); 
        sin_lambda = sin(lambda);

%  find necessary orbital elements of target at present (nu1) 
%  and future (nu2) locations 
        mu_m =  398600.4415; %*10^9;
        h_vec_tgt =  cross(r_tgt_RTN1,v_tgt_RTN1);
        p_tgt =  dot(h_vec_tgt,h_vec_tgt)/mu_m;
        ecc_vec_tgt =  cross(v_tgt_RTN1,h_vec_tgt)/mu_m-r_tgt_RTN1/mag_r_tgt;
        ecc_tgt =  norm(ecc_vec_tgt);
        a_tgt =  p_tgt/(1-ecc_tgt*ecc_tgt);
        % the check on ecc is sensitive. can result in diffs at certain
        % times
        if ecc_tgt > 0.0000001
          perigee_unit =  ecc_vec_tgt/ecc_tgt;
         else
          perigee_unit =  r_tgt_RTN1/mag_r_tgt;
        end;
        if ecc_tgt > 0.99
            r_tgt_RTN1
            v_tgt_RTN1
            r_int_RTN1
            v_int_RTN1
            dbstop;
        end    
        lambda_perigee =  atan2(perigee_unit(2,1),perigee_unit(1,1));
        nu1 =  -lambda_perigee;
        nu2 =  lambda-lambda_perigee;
        
%fprintf(1,'angles %11.6f %11.6f %11.6f %11.6f \n', phi*57.295, lambda*57.295, nu1*57.295, nu2*57.295);
        
%  find future position and velocity of target
        r2_tgt =  p_tgt / (1.0 + ecc_tgt*cos(nu2));
        P_vec =  perigee_unit;
        Q_vec =  cross([0 0 1]',P_vec);
        r2_vec_tgt =  r2_tgt*(cos(nu2)*P_vec + sin(nu2)*Q_vec);
        v2_vec_tgt =  sqrt(mu_m/p_tgt)*(-sin(nu2)*P_vec + (ecc_tgt+cos(nu2))*Q_vec);
        
%  rotate all to future target (RTN2) frame & adjust for phi
        rot_RTN1_to_RTN2 =  f_ECI_to_RTN_sal(r2_vec_tgt,v2_vec_tgt);
        r_tgt_RTN2 =  rot_RTN1_to_RTN2*r2_vec_tgt;
        v_tgt_RTN2 =  rot_RTN1_to_RTN2*v2_vec_tgt;
        
%  find interceptor SEZ components
        rot_to_SEZ =  zeros(3,3);
        rot_to_SEZ(1,1) = sin_phi*cos_lambda;
        rot_to_SEZ(1,2) = sin_phi*sin_lambda;
        rot_to_SEZ(1,3) = -cos_phi;
        rot_to_SEZ(2,1) = -sin_lambda;
        rot_to_SEZ(2,2) = cos_lambda;
        rot_to_SEZ(2,3) = 0.0;
        rot_to_SEZ(3,1) = cos_phi*cos_lambda;
        rot_to_SEZ(3,2) = cos_phi*sin_lambda;
        rot_to_SEZ(3,3) = sin_phi;
        r_int_SEZ = rot_to_SEZ*r_int_RTN1;
        v_int_SEZ = rot_to_SEZ*v_int_RTN1;

%  find position component positions  
        r_int_EQCM(1,1) = r_int_SEZ(3,1)-r_tgt_RTN2(1,1);

% try function for elliptic integral instead
        [ea0,m] = newtonnu ( ecc_tgt,nu1 );
        [ea1,m] = newtonnu ( ecc_tgt,nu2 );

        % fix quadrants for special cases
        if abs(ea1-ea0) > pi
            if ea0 < 0.0
                ea0 = 2.0*pi + ea0;
            else
                ea0 = 2.0*pi - ea0;
            end
        end        
        
        [F0, E0] = elliptic12( ea0, ecc_tgt^2 );
        [F1, E1] = elliptic12( ea1, ecc_tgt^2 );  % with implied tol of 1e-22, no error. Shows up at about 1e-5

        r_int_EQCM(2,1) = (a_tgt * (E1 - E0));  % arc length value
%        fprintf(1,'elliptic integral, %14.6f, %14.8f, %14.6f, %12.5f, %12.5f, ANS %12.6f  %12.6f \n', ...
%                      a_tgt, ecc_tgt, b, m0*rad, m1*rad, arclength, r_int_EQCM(2,1) );
   
        r_int_EQCM(3,1) = phi*r2_tgt;
%fprintf( 1,'%f,%f, ea0, %f, ea1, %f \n',E1-E0, ea1-ea0, ea0, ea1);    
%  find velocity component positions
        lamda_dot =  v_int_SEZ(2,1) / (mag_r_int*cos_phi);
        v_int_EQCM(1,1) = v_int_SEZ(3,1) - v_tgt_RTN2(1,1);
        v_int_EQCM(2,1) = lamda_dot*r2_tgt - (v_tgt_RTN1(2,1) / mag_r_tgt) * mag_r_tgt;
        v_int_EQCM(3,1) = (-v_int_SEZ(1,1) / mag_r_int)*r2_tgt;
        
        

      
      
      