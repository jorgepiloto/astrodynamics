% ------------------------------------------------------------------------------
%
%                           procedure biellip 
%
%  this procedure calculates the delta v's for a bi-elliptic transfer for either
%    circle to circle, or ellipse to ellipse.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    rinit       - initial position magnitude     er
%    r2          - interim orbit magnitude        er
%    rfinal      - final position magnitude       er
%    einit       - eccentricity of first orbit
%    efinal      - eccentricity of final orbit
%    nuinit      - true anomaly of first orbit    0 or pi rad
%    nufinal     - true anomaly of final orbit    0 or pi rad, opp of nuinit
%
%  outputs       :
%    deltava     - change in velocity at point a  er / tu
%    deltavb     - change in velocity at point b  er / tu
%    dttu        - time of flight for the trans   tu
%
%  locals        :
%    sme1        - mech energy of first orbit     er2 / tu
%    sme2        - mech energy of transfer orbit  er2 / tu
%    sme3        - mech energy of final orbit     er2 / tu
%    vinit       - velocity of first orbit at a   er / tu
%    vtransa     - velocity of trans orbit at a   er / tu
%    vtransb     - velocity of trans orbit at b   er / tu
%    vfinal      - velocity of final orbit at b   er / tu
%    ainit       - semimajor axis of first orbit  er
%    atrans      - semimajor axis of trans orbit  er
%    afinal      - semimajor axis of final orbit  er
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2007, 327, alg 37, ex 6-2
%function [deltava,deltavb,deltavc,dttu ] = biellip(rinit,rb,rfinal,einit,efinal,nuinit,nufinal);
% ----------------------------------------------------------------------------- }

function [deltava,deltavb,deltavc,dttu ] = biellip(rinit,rb,rfinal,einit,efinal,nuinit,nufinal);
     % --------------------  initialize values   ------------------- }
     mu = 1.0; % cannonical units

     ainit  = (rinit * (1.0 + einit * cos(nuinit))) / (1.0 - einit * einit );
     atran1 = (rinit + rb) * 0.5;
     atran2 = (rb + rfinal) * 0.5;
     afinal = (rfinal * (1.0 + efinal * cos(nufinal))) / (1.0 - efinal * efinal );

     deltava= 0.0;
     deltavb= 0.0;
     deltavc= 0.0;
     dttu   = 0.0;

     if ( einit < 1.0 ) & ( efinal < 1.0 ) 
 
     % -----------------  find delta v at point a  ----------------- }
         vinit  = sqrt( (2.0 * mu)/rinit - (mu/ainit) );
         vtran1a= sqrt( (2.0 * mu)/rinit - (mu/atran1) );
         deltava= abs( vtran1a - vinit );

     % -----------------  find delta v at point b  ----------------- }
         vtran1b= sqrt( (2.0 * mu)/rb - (mu/atran1) );
         vtran2b= sqrt( (2.0 * mu)/rb - (mu/atran2) );
         deltavb= abs( vtran1b - vtran2b );

     % -----------------  find delta v at point c  ----------------- }
         vtran2c= sqrt( (2.0 * mu)/rfinal - (mu/atran2) );
         vfinal = sqrt( (2.0 * mu)/rfinal - (mu/afinal) );
         deltavc= abs( vfinal - vtran2c );

     % ----------------  find transfer time of flight  ------------- }
         dttu= pi * sqrt( (atran1 * atran1 * atran1)/mu ) + ...
                pi * sqrt( (atran2 * atran2 * atran2)/mu );
     constastro;
     t1 = pi * sqrt( (atran1 * atran1 * atran1)/1 )*13.446852064;
     t2 = pi * sqrt( (atran2 * atran2 * atran2)/1 )*13.446852064;
     fprintf(1,' atran1   %11.7f  %11.7f km \n',atran1, atran1*re );
     fprintf(1,' atran2   %11.7f  %11.7f km \n',atran2, atran2*re );
     fprintf(1,' vinit   %11.7f  %11.7f km/s \n',vinit, vinit*velkmps );
     fprintf(1,' vtran1a  %11.7f  %11.7f km/s \n',vtran1a, vtran1a*velkmps );
     fprintf(1,' vtran1b  %11.7f  %11.7f km/s \n',vtran1b, vtran1b*velkmps );
     fprintf(1,' vtran2b  %11.7f  %11.7f km/s \n',vtran2b, vtran2b*velkmps );
     fprintf(1,' vtran2c  %11.7f  %11.7f km/s \n',vtran2c, vtran2c*velkmps );
     fprintf(1,' vfinal  %11.7f  %11.7f km/s \n',vfinal, vfinal*velkmps );
     fprintf(1,' t1  %11.7f t2 %11.7f min \n',t1, t2 );
     
     
     end;

