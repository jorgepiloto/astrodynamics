%
% ------------------------------------------------------------------------------
%
%                           function quadric
%
%  this function solves for the two roots of a quadric equation.  there are
%    no restrictions on the coefficients, and imaginary results are passed
%    out as separate values.  the general form is y = ax2 + bx + c.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  revisions
%    vallado     - convert to matlab              719-573-2600    3 dec 2002
%
%  inputs          description                    range / units
%    a           - coefficient of x squared term
%    b           - coefficient of x term
%    c           - constant
%    opt         - option for output              I all roots including imaginary
%                                                 R only real roots
%                                                 U only unique real roots (no repeated)
%
%  outputs       :
%    r1r         - real portion of root 1
%    r1i         - imaginary portion of root 1
%    r2r         - real portion of root 2
%    r2i         - imaginary portion of root 2
%
%  locals        :
%    discrim     - discriminate b2 - 4ac
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2007, 974
%
% [r1r,r1i,r2r,r2i] =  quadric   ( a,b,c,opt );
% ------------------------------------------------------------------------------  

function [r1r,r1i,r2r,r2i] =  quadric   ( a,b,c,opt );

        % --------------------  implementation   ----------------------
        small = 0.00000001;
        r1r = 0.0;
        r1i = 0.0;
        r2r = 0.0;
        r2i = 0.0;

        discrim = b*b - 4.0 *a*c;
%a
%b
%c
%discrim
        % ---------------------  real roots  --------------------------
        if ( abs(discrim) < small  )
            r1r = -b / ( 2.0 *a );
            r2r = r1r;
%            if (opt=='U')
%                r2r = 99999.9;
%              end;
          else
            if abs(a) < small
                 r1r = -c/b;
              else
            if ( discrim > 0.0  )
                r1r = ( -b + sqrt(discrim) ) / ( 2.0 *a );
                r2r = ( -b - sqrt(discrim) ) / ( 2.0 *a );
              else
                % ------------------ complex roots --------------------
                if (opt=='I')
                    r1r = -b / ( 2.0 *a );
                    r2r = r1r;
                    r1i =  sqrt(-discrim) / ( 2.0 *a );
                    r2i = -sqrt(-discrim) / ( 2.0 *a );
                  else
                    r1r = 99999.9;
                    r2r = 99999.9;
                  end;
              end;
            end;
          end;

