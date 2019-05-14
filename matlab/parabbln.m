% ------------------------------------------------------------------------------
%
%                           function parabbln
%
%  this function performs parabolic blending of an input zero crossing
%  function in order to find event times.
%
%  author        : david vallado                  719-573-2600    18 dec 2002
%
%  revisions
%                - fix eqt ref                                     3 jan 2003
%                - misc fixes                                      2 feb 2004
%
%  inputs          description                    range / units
%    p1,p2,p3    - function values used for blending
%
%  outputs       :
%    minfound    - test of success 
%    rootf       - root for the function
%    funrate     - function rate
%
%  locals        :
%
%  coupling      :
%    quadric     - find roots of a quadric
%
%  references    :
%    vallado       2007, 979
%
% [minfound,rootf,funrate] = parabbln( p1,p2,p3 );
% ------------------------------------------------------------------------------

function [minfound,rootf,funrate] = parabbln( p1,p2,p3 );

       rootf    = 0.0;
       funrate  = 0.0;
       minfound = 'n';
       
       % ------ set up function from C-37 --------
       aqd0 = p1;
       aqd1 = (-3.0*p1 + 4.0*p2 - p3)*0.5;
       aqd2 = (     p1 - 2.0*p2 + p3)*0.5;

       % --------------- solve roots of this function -------------
       opt = 'U';
       [r1r,r1i,r2r,r2i] = quadric( aqd2,aqd1,aqd0,opt );

       % ---------- search through roots to locate answers --------
       for indx2 = 1:2
           if (indx2==1)
               root = r1r;
             end;
           if (indx2==2)
               root = r2r;
             end;

           if (root >= 0.0) & (root <= 2.0)
%               [time]  = recovqd(t1,t2,t3, root); % should be 0.0!!!!!!
               [ans]  = recovqd(p1,p2,p3, root); % should be 0.0!!!!!!

               % ----- recover the function value derivative
               funrate = 2.0*aqd2*root + aqd1;

           end % if root between 0 and 2

         end % if indx2

