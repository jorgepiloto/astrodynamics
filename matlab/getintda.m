%
% -----------------------------------------------------------------------------
%
%                           function getintda
%
%  this function finds the integer equivalent of the 3 character string
%    representation of the day of the week.
%
%  author        : david vallado                  719-573-2600    5 jul 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    daystr      - day name string                'sun','mon' ...
%
%  outputs       :
%    dayn        - integer day equivalent         1 .. 7
%
%  locals        :
%    i           - index
%
%  coupling      :
%    none
%
% [dayn] = getintda( daystr );
% -----------------------------------------------------------------------------

function [dayn] = getintda( daystr );

        % ------------------------  implementation   --------------------------
        daytitle= str2mat('sun','mon','tue','wed','thr','fri','sat' );

        dayn = strmatch(daystr,daytitle);

