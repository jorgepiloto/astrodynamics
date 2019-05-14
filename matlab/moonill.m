% ------------------------------------------------------------------------------
%
%                           function moonill
%
%  this function calculates the illumination due to the moon.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    f           - moon hase angle                rad
%    moonel      - moon elevation                 rad
%
%  outputs       :
%    moonill     - moon illumination
%
%  locals        :
%                -
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2001, 295-297, eq 5-9
%
% [moonillum] = moonill ( f,moonel );
% ------------------------------------------------------------------------------

function [moonillum] = moonill ( f,moonel );

        % ------------------------  implementation   ------------------
        x= moonel/90.0;
        g= 1.0;

        if (moonel >= 20)
           l0= -1.95;
           l1=  4.06;
           l2= -4.24;
           l3=  1.56;
          elseif ((moonel >= 5.0) & (moonel < 20.0))
              l0=  -2.58;
             l1=  12.58;
             l2= -42.58;
             l3=  59.06;
           elseif ((moonel > -0.8) & (moonel < 5.0))
                 l0=   -2.79;
                 l1=   24.27;
                 l2= -252.95;
                 l3= 1321.29;
               else
                 l0= 0.0;
                 l1= 0.0;
                 l2= 0.0;
                 l3= 0.0;
                 f= 0.0;
                 g= 0.0;
               end

       l1= l0 + l1*x + l2*x*x + l3*x*x*x;
       l2= (-0.00868 *f - 2.2d-9*f*f*f*f);

%       hzparal =   0.9508 + 0.0518*cos( (134.9+477198.85*ttdb)*deg2rad )
%                + 0.0095*cos( (259.2-413335.38*ttdb)*deg2rad )
%                + 0.0078*cos( (235.7+890534.23*ttdb)*deg2rad )
%                + 0.0028*cos( (269.9+954397.70*ttdb)*deg2rad )   { deg }
%       hzparal  = realmod( hzparal*deg2rad, twopi )
%       l3= (2.0* power(10.0,(hzparal*rad / 0.951))*g ) { use g to eliminate neg el passes }

        moonillum= 10.0  ^ ( l1 + l2 );
        if ((moonillum < -1e+36) | (moonillum > 0.999 ))
            moonillum= 0.0;
          end

