% ------------------------------------------------------------------------------
%
%                           function rv2ntw
%
%  this function converts position and velocity vectors into normal (in-radial),
%    tangential (velocity), and normal (cross-track) coordinates. note that sometimes 
%    the first vector is called along-radial. the tangential direction is
%    always aligned with the velocity vector. this is the ntw system of
%    vallado. 
%
%  author        : david vallado                  719-573-2600    5 jul 2002
%
%  revisions
%                -
%  inputs          description                    range / units
%    r           - position vector                km
%    v           - velocity vector                km/s
%
%  outputs       :
%    rntw        - position vector                km
%    vntw        - velocity vector                km/s
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
% [rntw,vntw,transmat] = rv2ntw( r,v );
% ------------------------------------------------------------------------------

function [rntw, vntw, transmat] = rv2ntw(r, v);
        % compute satellite velocity vector magnitude
        vmag = mag(v);

        % in order to work correctly each of the components must be
        %  unit vectors !
        % in-velocity component
        tvec = v / vmag;

        % cross-track component
        wvec = cross(r, v);
        wvec = unit( wvec );

        % along-radial component
        nvec = cross(tvec, wvec);
        nvec = unit( nvec );

        % assemble transformation matrix from to ntw frame (individual
        %  components arranged in row vectors)
        transmat(1,1) = nvec(1);
        transmat(1,2) = nvec(2);
        transmat(1,3) = nvec(3);
        transmat(2,1) = tvec(1);
        transmat(2,2) = tvec(2);
        transmat(2,3) = tvec(3);
        transmat(3,1) = wvec(1);
        transmat(3,2) = wvec(2);
        transmat(3,3) = wvec(3);

        rntw = matvecmult(transmat,r,3);
        vntw = matvecmult(transmat,v,3);

