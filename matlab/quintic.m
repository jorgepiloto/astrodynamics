%
% ------------------------------------------------------------------------------
%
%                           function quintic
%
%  this function solves for the five roots of a quintic equation.  there are
%    no restrictions on the coefficients, and imaginary results are passed
%    out as separate values.  the general form is y  =  ax5 + bx4 + cx3 + dx2 +
%    ex + f.
%
%     this routine uses a newton-raphson search to find
%     the real roots (x) between 0 and 1 of the polynomial
%        a*x**5+b*x**4+...
%     the resulting 4th order equation is solved by calling root4,
%     if needed. (repeated roots are ignored)
%     taken from "numerical analysis" by maron, copyright 1982, pg 80
%
%  author        : david vallado                  719-573-2600   16 dec 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    a           - coefficient of x quintic term
%    b           - coefficient of x quartic term
%    c           - coefficient of x cubic term
%    d           - coefficient of x squared term
%    e           - coefficient of x term
%    f           - constant
%    opt         - option for output              I all roots including imaginary
%                                                 R only real roots
%                                                 U only unique real roots (no repeated)
%
%  outputs       :
%    r1r         - real portion of root 1
%    r1i         - imaginary portion of root 1
%    r2r         - real portion of root 2
%    r2i         - imaginary portion of root 2
%    r3r         - real portion of root 3
%    r3i         - imaginary portion of root 3
%    r4r         - real portion of root 4
%    r4i         - imaginary portion of root 4
%    r5r         - real portion of root 5
%    r5i         - imaginary portion of root 5
%
%  locals        :
%    temp1       - temporary value
%    temp2       - temporary value
%    p           - coefficient of x squared term where x cubed term is 1.0
%    q           - coefficient of x term where x cubed term is 1.0
%    r           - coefficient of constant term where x cubed term is 1.0
%    delta       - discriminator for use with cardans formula
%    e0          - angle holder for trigonometric solution
%    phi         - angle used in trigonometric solution
%    cosphi      - cosine of phi
%    sinphi      - sine of phi
%     f        polynomial coefficients
%     aa(5)       intermediate polynomial coefficients
%     dz          change in intermediate variable
%     f           value of polynomial function at root z
%     flag        flag for printing error message
%     fp          slope of polynomial function at root z
%     i           counter
%     maxdz       max step size for finding z
%     n           polynomial order
%     numrts      integer number of real roots
%     olddz       previous value of dz
%     small         tolerance
%     p5        real roots
%     xx(5)       intermediate real roots
%     z           intermediate root
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2001,
%
% [r1r,r1i,r2r,r2i,r3r,r3i,r4r,r4i,r5r,r5i]  =  quintic ( a,b,c,d,e,f,opt );
% ------------------------------------------------------------------------------

function [r1r,r1i,r2r,r2i,r3r,r3i,r4r,r4i,r5r,r5i]  =  quintic ( a,b,c,d,e,f,opt );

       % --------------------  implementation   ----------------------
       opt = 'U';
       rad        =  180.0/pi;
       onethird   =  1.0 /3.0;
       small      =  0.00000001;
       temp  =  1.0/24.0;
       r1r   =  0.0;
       r1i   =  0.0;
       r2r   =  0.0;
       r2i   =  0.0;
       r3r   =  0.0;
       r3i   =  0.0;

       % -----   test to see if there cannot be a crossing on the interval 0 to 1
       aa(5) = e;
       aa(4) = aa(5)+d;
       aa(3) = aa(4)+c;
       aa(2) = aa(3)+b;
       aa(1) = aa(2)+a;
       if ((f > 0.0) & (min(aa) > -f))
           return;
         end;
       if ((f < 0.0) & (max(aa) < -f))
           return;
         end;

       % ---------- check to see if it really is fifth order ----------
       if (abs(a) < small)
           [r1r,r1i,r2r,r2i,r3r,r3i,r4r,r4i] = quartic( a,b,c,d,e,opt );
           return
         end;
       if (abs(f) < small)
           [r1r,r1i,r2r,r2i,r3r,r3i,r4r,r4i] = quartic( b,c,d,e,f,opt );
           return
         end;

       % ----- find a good first guess for a real root between 0 and 1
       % ----- by assigning 5 evenly spaced points on the interval
       p1 = f;
       z = 0.25;
       p2 = ((((a*z+b)*z+c) *z+d) *z+e) *z + f;
       z = 0.5;
       p3 = ((((a*z+b)*z+c) *z+d) *z+e) *z + f;
       z = 0.75;
       p4 = ((((a*z+b)*z+c) *z+d) *z+e) *z + f;
       p5 = a+b+c+d+e+f;

       % -----   find polynomial coefficients (assumes the interval 0 to 4)
       aa =      p1;
       ab = (-50*p1 +96*p2 -72*p3 +32*p4  -6*p5)*temp;
       ac = ( 35*p1-104*p2+114*p3 -56*p4 +11*p5)*temp;
       ad = (-10*p1 +36*p2 -48*p3 +28*p4  -6*p5)*temp;
       ae = (    p1  -4*p2  +6*p3  -4*p4    +p5)*temp;

       % -----   find the roots of the fourth order polynomial and bound them
       [r1r,r1i,r2r,r2i,r3r,r3i,r4r,r4i] = quartic( aa,ab,ac,ad,ae,opt );
%       call bound4(xx,numrts,0.0,4.0);

       numrts = 0;
       if (abs(r1r) > small) & (abs(r1i) < small)
           numrts = numrts + 1;
         end;
       if (abs(r2r) > small) & (abs(r2i) < small)
           numrts = numrts + 1;
         end;
       if (abs(r3r) > small) & (abs(r3i) < small)
           numrts = numrts + 1;
         end;
       if (abs(r4r) > small) & (abs(r4i) < small)
           numrts = numrts + 1;
         end;

       % -----   if no real roots found   return to main routine
       if (numrts < 1)
           return
         end;
            
       % -----   initialize parameters using xx(1)/4 as a first guess for a real root
       % -----   (this rescales the root to the  interval 0 to 1)
       maxdz = 0.1;
       olddz = maxdz;
       z = r1r* 0.25;
       for i = 1,25
           f = ((((a*z+b)*z+c)*z+d)*z+e)*z+f;
           if (abs(f) < small)
               % ----- find roots of resulting 4th order equation and return
               a = a;
               b = b+z*aa(1);
               c = c+z*aa(2);
               d = d+z*aa(3);
               e = e+z*aa(4);
               [r1r,r1i,r2r,r2i,r3r,r3i,r4r,r4i] = quartic( a,b,c,d,e,opt );
               return
             end;
           % ----- find new dz and z
           fp = (((5*a*z+4*b)*z+3*c)*z+2*d)*z+e;
           if (abs(fp) > small)
               dz = -f/fp;
             end;
           if (abs(dz) > maxdz)
               dz = sign(maxdz) * dz;
               if (olddz*dz < 0.0)
                   maxdz = 0.8*maxdz;
                 end;
             end;
           olddz = dz;
           z = z+dz;
           % ----- test for z out of limits, if so there are no roots
           if ((z > 1.0) | (z <0.0))
               numrts = 0;
               return
             end;
         end; % ---- for

       % ---- once here, newton root search did not find an answer ----
       fprintf(1,' warning: no convergence in quintic ' );


