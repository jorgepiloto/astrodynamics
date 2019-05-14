% -----------------------------------------------------------------------------
%
%                           function hms2sec
%
%  this function converts hours, minutes and seconds into seconds from the
%    beginning of the day.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    hr          - hours                          0 .. 24
%    min         - minutes                        0 .. 59
%    sec         - seconds                        0.0 .. 59.99
%
%   outputs      :
%    utsec       - seconds                        0.0 .. 86400.0
%
%  locals        :
%    temp        - temporary variable
%
%  coupling      :
%    none.
%
% function [utsec ] = hms2sec( hr,min,sec );
% -----------------------------------------------------------------------------

function [utsec ] = hms2sec( hr,min,sec );

        % ------------------------  implementation   ------------------
        utsec  = hr * 3600.0 + min * 60.0 + sec;

