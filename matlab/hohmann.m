% ------------------------------------------------------------------------------
%
%                           procedure hohmann
%
%  this procedure calculates the delta v's for a hohmann transfer for either
%    circle to circle, or ellipse to ellipse.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    rinit       - initial position magnitude     er
%    rfinal      - final position magnitude       er
%    einit       - eccentricity of first orbit
%    efinal      - eccentricity of final orbit
%    nuinit      - true anomaly of first orbit    0 or pi rad
%    nufinal     - true anomaly of final orbit    0 or pi rad
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
%    vallado       2007, 327, alg 36, ex 6-1
%
%function [deltava,deltavb,dttu ] = hohmann (rinit,rfinal,einit,efinal,nuinit,nufinal);
% ----------------------------------------------------------------------------- }

function [deltava,deltavb,dttu ] = hohmann (rinit,rfinal,einit,efinal,nuinit,nufinal);
     % --------------------  initialize values   ------------------- }
     mu = 1.0; % cannonical units
     ainit  = (rinit * (1.0 + einit * cos(nuinit))) / (1.0 - einit * einit );
     atran  = ( rinit + rfinal ) / 2.0;
     afinal = (rfinal * (1.0 + efinal * cos(nufinal))) / (1.0 - efinal * efinal );
     deltava= 0.0;
     deltavb= 0.0;
     dttu   = 0.0;

     if ( einit < 1.0 ) || ( efinal < 1.0 )
         % -----------------  find delta v at point a  -------------- }
         vinit  = sqrt( (2.0 * mu)/rinit - (mu/ainit) );
         vtrana = sqrt( (2.0 * mu)/rinit - (mu/atran) );
         deltava= abs( vtrana - vinit );

         % -----------------  find delta v at point b  -------------- }
         vfinal = sqrt( (2.0 * mu)/rfinal - (mu/afinal) );
         vtranb = sqrt( (2.0 * mu)/rfinal - (mu/atran) );
         deltavb= abs( vfinal - vtranb );

         % ----------------  find transfer time of flight  ---------- }
         dttu= pi * sqrt( (atran * atran * atran) / mu );

     constastro;
     fprintf(1,' atran   %11.7f  %11.7f km \n',atran, atran*re );
     fprintf(1,' vinit   %11.7f  %11.7f km/s \n',vinit, vinit*velkmps );
     fprintf(1,' vtrana  %11.7f  %11.7f km/s \n',vtrana, vtrana*velkmps );
     fprintf(1,' vtranb  %11.7f  %11.7f km/s \n',vtranb, vtranb*velkmps );
     fprintf(1,' vfinal  %11.7f  %11.7f km/s \n',vfinal, vfinal*velkmps );
         
       end;

