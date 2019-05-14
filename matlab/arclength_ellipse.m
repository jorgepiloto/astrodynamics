function [arclength] = arclength_ellipse(a, b, theta0, theta1)
%ARCLENGTH_ELLIPSE Calculates the arclength of ellipse.
%
%   ARCLENGTH_ELLIPSE(A, B, THETA0, THETA1) Calculates the arclength of ellipse 
%   using the precise formulas based on the representation of 
%   the arclength by the Elliptic integral of the second kind.
%
%   Ellipse parameters:
%       T - measured in radians from 0 in the positive direction, 
%           Period: 2*Pi
%       A - major axis
%       B - minor axis
%   
%   Parametric equations:
%       x(t) = a.cos(t)
%       y(t) = b.sin(t)
%
%   Cartesian equation:
%   x^2/a^2 + y^2/b^2 = 1
%
%   Eccentricity:
%       e = Sqrt(1 - (a/b)^2)
%
%   Focal parameter:
%       b^2/Sqrt(a^2 - b^2)
%
%   Foci:
%       (-Sqrt(a^2 - b^2), 0)   OR   (Sqrt(a^2 - b^2), 0)
%
%   Arclength:
%       b*EllipticE( t, 1 - (a/b)^2 )
%
%   Mathematica Test 1:
%       In:= b = 10; a = 5;
%            SetPrecision[b*EllipticE[2Pi, 1.0- a^2/b^2],20]
%      Out:= 48.442241102738385905
%
%   Mathematica Test 2:
%       In:= b = 10; a = 5;
%            SetPrecision[b*(EllipticE[Pi/2-Pi/10, 1.0- a^2/b^2]-EllipticE[Pi/10, 1.0- a^2/b^2]),20]
%      Out:= 7.3635807913930495516
%
%   MATLAB Test 1:
%       % full ellipse
%       arclength = arclength_ellipse(5,10)
%       arclength =
%           48.442241102738436
%
%   MATLAB Test 2:
%       % arclength ellipse
%       arclength = arclength_ellipse(5,10,pi/10,pi/2)
%       arclength =
%           7.363580791393055
%
%   References:
%   @see http://mathworld.wolfram.com/Ellipse.html
%   @see http://www.wolframalpha.com/input/?i=ellipse+arc+length&lk=1&a=ClashPrefs_*PlaneCurve.Ellipse.PlaneCurveProperty.ArcLength-
%

% Special thanks to for bug correction
%    drbitboy (Brian Carcich) https://github.com/drbitboy
% 2015-07-14 (New Horizons flyby of Pluto)
%
% 1) Old code returned values that were in error
% 1.1)  arclength_ellipse(1., .5, pi*.001, pi*.002) returned 0
% 1.2)  arclength_ellipse(1., .5, pi*.002, pi*.001) returned -.0003*pi instead of pi correct .0005*pi
% 1.3)  arclength_ellipse(1., .5, theta0, theta1) did not return the negative of the same call with the thetas reversed
% 2) Angles theta0 and theta1 were always interpreted as measured from the semi-minor axis
%
% 3) Corrected code:
% 3.1) Angle theta is measured from the positive a axis
% 3.2) The standard form of the b*E(phi,m) arc length integral has m = 1 - (a/b)^2
% 3.2.1) N.B. That only only works if b>a
% 3.3) If a>b, then an alternate formula is used:  a*E(PI/2 - phi, m') where m' = 1 - (b/a)^2
% 3.4) A few simple cases will show that the new code is correct
%        arclength_ellipse(1, .5, pi*.001, pi*.002) ~  pi*.0005
%        arclength_ellipse(1, .5, pi*.002, pi*.001) = -arclength(1, .5, pi*.001, pi*.002) ~ -pi*.0005
%        arclength_ellipse(1., 2., pi*.001, pi*.002) ~ pi*.002
%        arclength_ellipse(1, .5, pi/2 - pi*.002, pi/2 - pi*.001) ~ -pi*.001
%        arclength_ellipse(1, 2., pi/2 - pi*.002, pi/2 - pi*.001) ~ -pi*.001
%        etc.

% Copyright Elliptic Project 2011
% For support,  
%     moiseev.igor[at]gmail.com
%     Moiseev Igor 

%arguments
if nargin ~= 2 && nargin ~= 4,
 error('ARCLENGTH_ELLIPSE: Requires two or four inputs.')
 return
end

if nargin == 2,
 theta0 = 0;
 theta1 = 2*pi;
end

% Default solution for a==b (circles)
arclength = a.*(theta1-theta0);

% Ellipses (a<b or a>b)
if(a<b)
    % Theta measured from a axis = semi-MINOR axis
    % Use standard formulation for E(phi,m)
    [F1, E1] = elliptic12( theta1, 1 - (a./b).^2 );
    [F0, E0] = elliptic12( theta0, 1 - (a./b).^2 );
    arclength = b.*(E1 - E0);
elseif(a>b)   
    % Theta measured from a axis = semi-MAJOR axis
    % Standard formulation will not work ((1-(a/b)^2) < 0); instead use PI/2 - phi and b/a instead of a/b
    [F1, E1] = elliptic12( pi/2 - theta1, 1 - (b./a).^2 );
    [F0, E0] = elliptic12( pi/2 - theta0, 1 - (b./a).^2 );
    % d(PI/2 - phi)/dphi = -1, so reverse operands in this difference to flip sign:
    arclength = a.*(E0 - E1);
end

return;


