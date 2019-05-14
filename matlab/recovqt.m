%
% ------------------------------------------------------------------------------
%
%                           function recovqt
%
%  this function recovers the time and function values in quartic blending routines.
%
%  author        : david vallado                  719-573-2600    11 dec 2002
%
%  revisions
%                - misc fixes                                      2 feb 2004
%
%  inputs          description                    range / units
%    p1,p2,p3,p4,p5,p6
%                - function values used for blending
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
%    vallado       2001, 900
%
% [funvalue] = recovqt ( p1,p2,p3,p4,p5,p6, root );
% ------------------------------------------------------------------------------

function [funvalue] = recovqt ( p1,p2,p3,p4,p5,p6, root );

       % ------ set up function from C-45 --------
       %  aqit5*x**5 + aqit4*x**4 + etc
       temp = 1.0/24.0;
       aqit0 =                   p3;
       aqit1 = ( 2*p1 -16*p2        +16*p4  -2*p5)*temp;
       aqit2 = (-1*p1 +16*p2 -30*p3 +16*p4    -p5)*temp;
       aqit3 = (-9*p1 +39*p2 -70*p3 +66*p4 -33*p5  +7*p6)*temp;
       aqit4 = (13*p1 -64*p2+126*p3-124*p4 +61*p5 -12*p6)*temp;
       aqit5 = (-5*p1 +25*p2 -50*p3 +50*p4 -25*p5  +5*p6)*temp;

       % ----- recover the variable value
       funvalue = aqit5*root^5 +aqit4*root^4 + aqit3*root^3 +...
                  aqit2*root^2 + aqit1*root + aqit0;


