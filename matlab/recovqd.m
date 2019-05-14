%
% ------------------------------------------------------------------------------
%
%                           function recovqd
%
%  this function recovers the time and function values in parabolic blending routines.
%
%  author        : david vallado                  719-573-2600     3 jan 2003
%
%  revisions
%                - misc fixes                                      2 feb 2004
%
%  inputs          description                    range / units
%    p1,p2,p3    - function values used for blending
%    root        - root used as variable
%
%  outputs       :
%    funvalue    - function value
%
%  locals        :
%    none
%
%  coupling      :
%    none
%
%  references    :
%    vallado       2001, 895
%
% [funvalue] = recovqd ( p1,p2,p3, root );
% ------------------------------------------------------------------------------

function [funvalue] = recovqd ( p1,p2,p3, root );

       % ------ set up function from C-37 --------
       aqd0 = p1;
       aqd1 = (-3.0*p1 + 4.0*p2 - p3)*0.5;
       aqd2 = (     p1 - 2.0*p2 + p3)*0.5;

       % ----- recover the variable value
       funvalue = aqd2*root^2 + aqd1*root + aqd0;



