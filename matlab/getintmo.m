%
% -----------------------------------------------------------------------------
%
%                           function getintmon
%
%  this function finds the integer equivalent of the 3 character string
%    representation of month.
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    monstr      - month name                     'jan','feb' ...
%
%  outputs       :
%    mon         - integer month equivalent       1 .. 12
%
%  locals        :
%    i           - index
%
%  coupling      :
%    none
%
% [mon] = getintmon( monstr );
% -----------------------------------------------------------------------------

function [mon] = getintmon( monstr );

        % ------------------------  implementation   --------------------------
        monthtitle= str2mat('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec');

        mon = strmatch(monstr,monthtitle);


