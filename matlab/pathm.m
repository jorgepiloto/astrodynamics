% ------------------------------------------------------------------------------
%
%                           function pathm
%
%  this function determines the end position for a given range and azimuth
%    from a given point.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    llat        - start geocentric latitude      -pi/2 to  pi/2 rad
%    llon        - start longitude (west -)       0.0  to 2pi rad
%    range       - range between points           er
%    az          - azimuth                        0.0  to 2pi rad
%
%  outputs       :
%    tlat        - end geocentric latitude        -pi/2 to  pi/2 rad
%    tlon        - end longitude (west -)         0.0  to 2pi rad
%
%  locals        :
%    sindeltan   - sine of delta n                rad
%    cosdeltan   - cosine of delta n              rad
%    deltan      - angle between the two points   rad
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2001, 774-776, eq 11-6, eq 11-7
%
% [tlat,tlon] = pathm ( llat, llon, range, az );
% ------------------------------------------------------------------------------

function [tlat,tlon] = pathm ( llat, llon, range, az );

        twopi =     2.0*pi;

        % -------------------------  implementation   -----------------
        small =     0.00000001;

        az= rem( az,twopi );
        if ( llon < 0.0  )
            llon= twopi + llon;
          end
        if ( range > twopi )
            range= rem( range,twopi );
          end

        % ----------------- find geocentric latitude  -----------------
        tlat = asin( sin(llat)*cos(range) + cos(llat)*sin(range)*cos(az) );

        % ---- find delta n, the angle between the points -------------
        if ( (abs(cos(tlat)) > small) & (abs(cos(llat)) > small) )
            sindn = sin(az)*sin(range) / cos(tlat);
            cosdn = ( cos(range)-sin(tlat)*sin(llat) ) / ( cos(tlat)*cos(llat) );
            deltan= atan2(sindn,cosdn);
          else
            % ------ case where launch is within 3nm of a pole --------
            if ( abs(cos(llat)) <= small )
                if ( (range > pi) & (range < twopi) )
                    deltan= az + pi;
                  else
                    deltan= az;
                  end
              end
            % ----- case where end point is within 3nm of a pole ------
            if ( abs( cos(tlat) ) <= small )
                deltan= 0.0;
              end
          end

        tlon= llon + deltan;
        if ( abs(tlon) > twopi )
            tlon= rem( tlon,twopi );
          end
        if ( tlon < 0.0  )
            tlon= twopi + tlon;
          end

