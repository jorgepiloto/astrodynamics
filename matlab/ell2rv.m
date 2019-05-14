%
% ecliptic latitude longitude to position and velocity
% dav 28 mar 04
%
% [rijk,vijk] = ell2rv ( rr,ecllon,ecllat,drr,decllon,decllat );

function [rijk,vijk] = ell2rv ( rr,ecllon,ecllat,drr,decllon,decllat );

        % --------------------  implementation   ----------------------
        obliquity= 0.40909280;   %23.439291 /rad

        r(1)= rr*cos(ecllat)*cos(ecllon);
        r(2)= rr*cos(ecllat)*sin(ecllon);
        r(3)= rr*sin(ecllat);

        v(1)= drr*cos(ecllat)*cos(ecllon) ...
                 - rr*sin(ecllat)*cos(ecllon)*decllat ...
                 - rr*cos(ecllat)*sin(ecllon)*decllon;
        v(2)= drr*cos(ecllat)*sin(ecllon) ...
                 - rr*sin(ecllat)*sin(ecllon)*decllat ...
                 + rr*cos(ecllat)*cos(ecllon)*decllon;
        v(3)= drr*sin(ecllat) + rr*cos(ecllat)*decllat;

        [rijk] = rot1 ( r, -obliquity );
        [vijk] = rot1 ( v, -obliquity );



