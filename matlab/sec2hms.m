% -----------------------------------------------------------------------------
%
%                           function sec2hms
%
%  this function converts seconds from the beginning of the day into hours,
%    minutes and seconds.
%
%  author        : david vallado                  719-573-2600   25 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    utsec       - seconds                        0.0 .. 86400.0
%
%  outputs       :
%    hr          - hours                          0 .. 24
%    min         - minutes                        0 .. 59
%    sec         - seconds                        0.0 .. 59.99
%
%  locals        :
%    temp        - temporary variable
%
%  coupling      :
%    none.
%
% [hr,min,sec] = sec2hms( utsec );
% -----------------------------------------------------------------------------

function [hr,min,sec] = sec2hms( utsec );

        % ------------------------  implementation   ------------------
        temp  = utsec / 3600.0;
        hr    = fix( temp );
        min   = fix( (temp - hr)* 60.0 );
        sec   = (temp - hr - min/60.0 ) * 3600.0;

