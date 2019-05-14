function invE = inverselliptic2(E,m,tol)
% INVERSELLIPTIC2 evaluates the value of the INVERSE Incomplete Elliptic Integrals 
% of the Second Kind.
%
% INVE = INVERSELLIPTIC2(E,M,TOL) where E is a value of the integral to 
% be inverted, 0<M<1 is the module and TOL is the tolerance (optional). 
% Default value for the tolerance is eps = 2.220e-16.
%
% INVERSELLIPTIC2 uses the method described by Boyd J. P. 
% to determine the value of the inverse Incomplete Elliptic Integrals 
% of the Second Kind using the “Empirical” initialization to 
% the Newton’s iteration method [1]. 
%
% NOTICE. Please pay attention to the definition of the elliptic functions
% which follows the Abramovitz et al [2], for more theory on elliptic
% functions please consult the Lawden book [3].
%
% Elliptic integral of the second kind:
%
% E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
%
% “Empirical” initialization [1]:
%
% T0(z,m) = pi/2 + sqrt(r)/(theta ? pi/2)
%
% where 
% z \in [?E(pi/2,m), E(pi/2,m)]x[0, 1], value of the entire parameter space
% r = sqrt((1-m)^2 + zeta^2)
% zeta = 1 - z/E(pi/2,m)
% theta = atan((1 - m)/zeta)
%
%
% Example:
% % modulus and phase in degrees
% [phi,alpha] = meshgrid(0:5:90, 0:2:90);
% % values of integrals
% [F,E] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);
% % values of inverse 
% invE = inverselliptic2(E, sin(pi/180*alpha).^2);
% % the difference between phase phi and invE should close to zero
% phi - invE * 180/pi
%
% See also ELLIPKE, ELLIPTIC12.
%
% References:
% [1] J. P. Boyd, "Numerical, Perturbative and Chebyshev Inversion 
% of the Incomplete Elliptic Integral of the Second Kind", Applied Mathematics and Computation (January 2012)
% [2] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions", 
% Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
% [3] D. F. Lawden, "Elliptic Functions and Applications"
% Springer-Verlag, vol. 80, 1989
%
% Copyright (C) 2011 by Elliptic Project. All rights reserved.

% GNU GENERAL PUBLIC LICENSE Version 2, June 1991
% http://www.gnu.org/licenses/gpl.html 
% Everyone is permitted to copy and distribute verbatim copies of this 
% script under terms and conditions of GNU GENERAL PUBLIC LICENSE. 
% 
% For support, please reply to 
% moiseev.igor[at]gmail.com
% Moiseev Igor
%
% ELLIPTIC PROJECT: http://elliptic.googlecode.com
% Group: 

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(E) || ~isreal(m)
error('Input arguments must be real.');
end

if length(m)==1, m = m(ones(size(E))); end
if length(E)==1, E = E(ones(size(m))); end
if ~isequal(size(m),size(E)), error('E and M must be the same size.'); end

invE = zeros(size(E)); 

% make a row vector
m = m(:);
E = E(:);

if any(m < 0) || any(m > 1), error('M must be in the range 0 <= M <= 1.'); end
% cdav change for small eccentricities
if abs(m) < 1e-7 
    m=1e-7;
end 
% inputs
z = E; mu = 1-m;

% complete integral initialization
[~,E1] = ellipke(m,tol); 

zeta = 1 - z./E1;
r = sqrt(zeta.*zeta+mu.*mu); 
theta = atan(mu./(z+eps));

% “Empirical” initialization [1]
invE(:) = pi/2 + sqrt(r).*(theta - (pi/2)); 

for iter=1:4
  [~, E] = elliptic12(invE(:),m,tol);
  invE(:) = invE(:)-(E - z)./sqrt( 1-m.*sin(invE(:)).^2 ); 
end
return;


