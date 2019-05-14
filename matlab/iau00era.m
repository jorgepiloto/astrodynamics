%
% ----------------------------------------------------------------------------
%
%                           function iau00era
%
%  this function calulates the transformation matrix that accounts for the
%    effects of sidereal time via the earth rotation angle.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  inputs          description                    range / units
%    jdut1       - julian date of ut1             days
%
%  outputs       :
%    st          - transformation matrix for pef-ire
%
%  locals        :
%    tdut1       - julian centuries of ut1        days
%    era         - earth rotation angle           rad
%
%  coupling      :
%
%  references    :
%    vallado       2004, 212
%
% [st]  = iau00era (jdut1 );
% ----------------------------------------------------------------------------

function [st]  = iau00era (jdut1 );

        sethelp;

        constastro;

        % julian centuries of ut1
        tut1d= jdut1 - 2451545.0;

        era = twopi * ( 0.7790572732640 + 1.00273781191135448 * tut1d );
        era = rem (era,twopi);

        if iauhelp == 'y'
            fprintf(1,'era%11.7f  \n',era*180/pi );
          end;

        % transformation matrix
        st(1,1) =  cos(era);
        st(1,2) = -sin(era);
        st(1,3) =  0.0;

        st(2,1) =  sin(era);
        st(2,2) =  cos(era);
        st(2,3) =  0.0;

        st(3,1) =  0.0;
        st(3,2) =  0.0;
        st(3,3) =  1.0;


