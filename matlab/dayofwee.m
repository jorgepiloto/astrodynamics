% -----------------------------------------------------------------------------
%
%                           function dayofwee
%
%  this function finds the day of the week. integers are used for the days,
%    1 = 'sun', 2 = 'mon', ... 7 = 'sat'.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    jd          - julian date of interest        days from 4713 bc
%
%  outputs       :
%    dayofweek   - answer                         1 to 7
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2007, 188, eq 3-39
%
% dayofweek = dayofwee(jd);
% -----------------------------------------------------------------------------

function dayofweek = dayofwee(jd);

        % ------------------------  implementation   ------------------
        % ------- be sure jd is at 0.0d0 h on the day of interest -----
        jd = floor(jd + 0.5);

        dayofweek = fix( jd - 7 * fix( (jd+1)/7 ) + 2 );

