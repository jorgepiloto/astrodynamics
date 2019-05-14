% ------------------------------------------------------------------------------
%
%                           function rv2rsw
%
%  this function converts position and velocity vectors into radial, tangential (in-
%    track), and normal (cross-track) coordinates. note that there are numerous 
%    nomenclatures for these systems. this is the rsw system of vallado. the reverse
%    values are found using the transmat transpose. 
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%                -
%  inputs          description                    range / units
%    reci        - position vector                km
%    veci        - velocity vector                km/s
%
%  outputs       :
%    rrsw        - position vector                km
%    vrsw        - velocity vector                km/s
%
%  locals        :
%    temp        - temporary position vector
%
%  coupling      :
%    mag         - magnitude of a vector
%    matvecmult  - multiply a square matrix by a vector (single column)
%
%  references    :
%    vallado       2007, 172
%
% [rrsw,vrsw,transmat] = rv2rsw( reci,veci );
% ------------------------------------------------------------------------------

function [rrsw, vrsw, transmat] = rv2rsw(reci, veci);

        % each of the components must be unit vectors
        % radial component
        rvec = unit(reci);

        % cross-track component
        wvec    = cross(reci, veci);
        wvec    = unit( wvec );

        % along-track component
        svec    = cross(wvec, rvec);
        svec    = unit( svec );

        % assemble transformation matrix from eci to rsw frame (individual
        % components arranged in row vectors)
        transmat(1,1) = rvec(1);
        transmat(1,2) = rvec(2);
        transmat(1,3) = rvec(3);
        transmat(2,1) = svec(1);
        transmat(2,2) = svec(2);
        transmat(2,3) = svec(3);
        transmat(3,1) = wvec(1);
        transmat(3,2) = wvec(2);
        transmat(3,3) = wvec(3);

        rrsw = matvecmult(transmat, reci, 3);
        vrsw = matvecmult(transmat, veci, 3);

%   alt approach
%       rrsw(1) = mag(reci);
%       rrsw(2) = 0.0;
%       rrsw(3) = 0.0;
%       vrsw(1) = dot(reci,veci)/rrsw(1);
%       vrsw(2) = sqrt(veci(1)^2 + veci(2)^2 + veci(3)^2 - vrsw(1)^2);
%       vrsw(3) = 0.0;

