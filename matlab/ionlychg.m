% ------------------------------------------------------------------------------
%
%                           procedure ionlychg
%
%  this procedure calculates the delta v's for a change in inclination only.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    deltai      - change in inclination          rad
%    vinit       - initial velocity vector        er/tu
%    fpa         - flight path angle              rad
%
%  outputs       :
%    deltavionly - answer
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2007, 346, alg 39, ex 6-4
%function [ deltavionly] = ionlychg(deltai,vinit,fpa);
% ----------------------------------------------------------------------------- }

function [ deltavionly] = ionlychg(deltai,vinit,fpa);

       deltavionly = 2.0 * vinit * cos(fpa) * sin(0.5 * deltai);

