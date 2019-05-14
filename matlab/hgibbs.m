%
% ------------------------------------------------------------------------------
%
%                           function hgibbs
%
%  this function implements the herrick-gibbs approximation for orbit
%    determination, and finds the middle velocity vector for the 3 given
%    position vectors.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    r1          - ijk position vector #1         km
%    r2          - ijk position vector #2         km
%    r3          - ijk position vector #3         km
%    use seconds to provide more accuracy
%    t1         - time julian date of 1st sighting    days from 4713 bc sec
%    t2         - time julian date of 2nd sighting    days from 4713 bc sec
%    t3         - time julian date of 3rd sighting    days from 4713 bc sec
%
%  outputs       :
%    v2          - ijk velocity vector for r2     km / s
%    theta       - angl between vectors          rad
%    error       - flag indicating success        'ok',...
%
%  locals        :
%    dt21        - time delta between r1 and r2   s
%    dt31        - time delta between r3 and r1   s
%    dt32        - time delta between r3 and r2   s
%    p           - p vector    r2 x r3
%    pn          - p unit vector
%    r1n         - r1 unit vector
%    theta1      - temporary angl between vec    rad
%    tolangle    - tolerance angl  (1 deg)       rad
%    term1       - 1st term for hgibbs expansion
%    term2       - 2nd term for hgibbs expansion
%    term3       - 3rd term for hgibbs expansion
%    i           - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    cross       - cross product of two vectors
%    dot         - dot product of two vectors
%    unit        - unit vector
%    lncom3      - combination of three scalars and three vectors
%    angl       - angl between two vectors
%
%  references    :
%    vallado       2007, 462, alg 52, ex 7-4
%
% [v2, theta,theta1,copa, error ] =  hgibbs ( r1,r2,r3,jd1,jd2,jd3 );
% ------------------------------------------------------------------------------

function [v2, theta,theta1,copa, error ] =  hgibbs ( r1,r2,r3,jd1,jd2,jd3 );

% -------------------------  implementation   -------------------------
        constmath;
        constastro;

        error =  '          ok';
        theta = 0.0;
        theta1= 0.0;
        magr1 = mag( r1 );
        magr2 = mag( r2 );
        magr3 = mag( r3 );
        for i= 1 : 3
            v2(i)= 0.0;
          end

        tolangle= 0.01745329251994;
        dt21= (jd2-jd1)*86400.0;
        dt31= (jd3-jd1)*86400.0;    % differences in times in secs
        dt32= (jd3-jd2)*86400.0;

        p = cross( r2,r3 );
        pn = unit( p );
        r1n = unit( r1 );
        copa=  asin( dot( pn,r1n ) );
        if ( abs( dot(r1n,pn) ) > 0.017452406 )
            error= 'not coplanar';
          end

        % --------------------------------------------------------------
%       check the size of the angles between the three position vectors.
%       herrick gibbs only gives "reasonable" answers when the
%       position vectors are reasonably close.  10 deg is only an estimate.
        % --------------------------------------------------------------
        theta  = angl( r1,r2 );
        theta1 = angl( r2,r3 );
        if ( (theta > tolangle) | (theta1 > tolangle) )  
            error= '   angl > 1ø';
          end

        % ----------- perform herrick-gibbs method to find v2 ---------
        term1= -dt32*( 1.0/(dt21*dt31) + mu/(12.0*magr1*magr1*magr1) );
        term2= (dt32-dt21)*( 1.0/(dt21*dt32) + mu/(12.0*magr2*magr2*magr2) );
        term3=  dt21*( 1.0/(dt32*dt31) + mu/(12.0*magr3*magr3*magr3) );

        v2 =  term1*r1 + term2* r2 + term3* r3;

%     fprintf(1,'p     %11.7f   %11.7f  %11.7f km2 \n',p);
%     fprintf(1,'theta     %11.7f   %11.7f deg \n',theta*180/pi, theta1*180/pi );

