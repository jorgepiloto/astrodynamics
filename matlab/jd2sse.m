% -----------------------------------------------------------------------------
%
%                           function jd2sse.m
%
%  this function finds the seconds since epoch (1 Jan 2000) given the julian date
%
%  author        : david vallado                  719-573-2600   12 dec 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    jd          - julian date                    days from 4713 bc
%
%  outputs       :
%    sse         - seconds since epoch 1 jan 2000
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
%  references    :
%    none.
%
% sse = jd2sse( jd );
% -----------------------------------------------------------------------------

function sse = jd2sse( jd );

        % ------------------------  implementation   ------------------
        sse = (jd - 2451544.5) * 86400.0;


