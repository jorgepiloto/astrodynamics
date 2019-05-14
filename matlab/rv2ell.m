%
% position and velocity to ecliptic latitude longitude
% dav 28 mar 04
%
% [rr,ecllon,ecllat,drr,decllon,decllat] = rv2ell (rijk,vijk);

function [rr,ecllon,ecllat,drr,decllon,decllat] = rv2ell (rijk,vijk);

        % --------------------  implementation   ----------------------
        constmath;
        obliquity= 0.40909280;   %23.439291 /rad

        [r] = rot1 ( rijk, obliquity );
        [v] = rot1 ( vijk, obliquity );

        % ------------- calculate angles and rates ----------------
        rr= mag(r);
        temp= sqrt( r(1)*r(1) + r(2)*r(2) );
        if ( temp < small )
            temp1= sqrt( v(1)*v(1) + v(2)*v(2) );
            if ( abs(temp1) > small )
                ecllon= atan2( v(2) , v(1) );
              else
                ecllon= 0.0;
              end;
          else
            ecllon= atan2( r(2) , r(1) );
          end;
        ecllat= asin( r(3)/rr );

        temp1= -r(2)*r(2) - r(1)*r(1);  % different now
        drr= dot(r,v)/rr;
        if ( abs( temp1 ) > small )
            decllon= ( v(1)*r(2) - v(2)*r(1) ) / temp1;
          else
            decllon= 0.0;
          end;
        if ( abs( temp ) > small )
            decllat= ( v(3) - drr*sin( ecllat ) ) / temp;
          else
            decllat= 0.0;
          end;

