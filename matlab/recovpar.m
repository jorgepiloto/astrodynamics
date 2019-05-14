%
% ------------------------------------------------------------------------------
%
%                           function recovpar
%
%  this function recovers the time and function values in parabolic blending routines.
%
%  author        : david vallado                  719-573-2600    11 dec 2002
%
%  revisions
%                - fix for single parabbln                        19 sep 2003
%                - misc fixes                                      2 feb 2004
%
%  inputs          description                    range / units
%    p1,p2,p3,p4 - function values used for blending
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
%    vallado       2001, 897
%
% [funvalue] = recovpar ( p1,p2,p3,p4, root );
% ------------------------------------------------------------------------------

function [funvalue] = recovpar ( p1,p2,p3,p4, root );

       % ------ set up function from C-39 -------
       %  acut3*x**3 + acut2*x**2 + etc
       acut0 = p2;
       acut1 = (-p1 + p3)*0.5;
       acut2 =      p1 - 2.5*p2 + 2.0*p3 - 0.5*p4;
       acut3 = -0.5*p1 + 1.5*p2 - 1.5*p3 + 0.5*p4;

       % ----- recover the variable value
       funvalue = acut3*root^3 + acut2*root^2 + acut1*root + acut0;


