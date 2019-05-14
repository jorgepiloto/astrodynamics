% ------------------------------------------------------------------------------
%
%                           function f_ECI_to_RTN_sal
%
%  this function converts finds rotation matrix from ECI to RTN
%
%  author        : sal alfano             719-573-2600    11 aug 2010
% ------------------------------------------------------------------------------

function rot_ECI_to_RTN = f_ECI_to_RTN_sal( r_ECI,v_ECI );

        R_unit=r_ECI/norm(r_ECI);
        h=cross(r_ECI,v_ECI);
        N_unit=h/norm(h);
        T_unit=cross(N_unit,R_unit);
        rot_ECI_to_RTN=horzcat(R_unit,T_unit,N_unit)';
