% ------------------------------------------------------------------------------
%
%                           procedure onetang
%
%  this procedure calculates the delta v's for a one tangent transfer for either
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
%    nu2         - true anomaly of second orbit   same quad as nuinit, rad
%    nufinal     - true anomaly of final orbit    0 or pi rad
%
%  outputs       :
%    deltava     - change in velocity at point a  er / tu
%    deltavb     - change in velocity at point b  er / tu
%    dttu        - time of flight for the transf  tu
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
%    e           - ecc anomaly of trans at b      rad
%    ratio       - ratio of initial to final
%                    orbit radii
%
%  coupling      :
%    atan2       - arc tangent rountine that solves quadrant ambiguities
%
%  references    :
%    vallado       2007, 335, alg 38, ex 6-3
%function [deltava, deltavb, dttu, etran, atran, vtrana, vtranb ] = onetang(rinit, rfinal, einit, efinal, nuinit, nutran);
% ----------------------------------------------------------------------------- }

function [deltava, deltavb, dttu, etran, atran, vtrana, vtranb ] = onetang(rinit, rfinal, einit, efinal, nuinit, nutran);
     % --------------------  initialize values   ------------------- }
     mu = 1.0; % cannonical units
     re = 6378.137;

     e = 0.0;
     deltava= 0.0;
     deltavb= 0.0;
     dttu    = 0.0;
     ratio  = rinit/rfinal;
     if abs(nuinit) < 0.01 % check 0 or 180 }
         etran  = (( ratio - 1.0 ) / ( cos(nutran) - ratio ));  % init at perigee 
         etran  = (( -rfinal + rinit ) / ( rfinal * cos(nutran) - rinit ));  % init at perigee 
         eainit= 0.0;
       else
         etran  = (( ratio - 1.0 ) / ( cos(nutran) + ratio ));  % init at perigee 
         etran  = (( -rfinal + rinit ) / ( rfinal * cos(nutran) + rinit ));  % init at apogee 
         eainit= pi;
       end;
     if etran >= 0.0  
         ainit = (rinit * (1.0 + einit * cos(nuinit))) / (1.0 - einit * einit );
         afinal= (rfinal * (1.0 + efinal * cos(nutran))) / (1.0 - efinal * efinal );
                      % nutran is used since it = nufinal!! }
%    fprintf(1,' ainti and final   %11.7f  %11.7f km \n',ainit*re, afinal*re );
%ainit = rinit;
%afinal= rfinal;
         if abs( etran-1.0 ) > 0.000001  
             if abs(nuinit) < 0.01 % check 0 or 180 }
                 atran = (rinit * (1.0 + etran * cos(nuinit))) / (1.0 - etran*etran ); % per }
               else
                 atran = (rinit * (1.0 + etran * cos(nuinit))) / (1.0 + etran * etran ); %  apo }
                 atran= rinit/(1.0 + etran);
               end
           else
             atran = 999999.9;  % infinite for parabolic orbit }
           end;

         ptran = rinit * ( 1.0 + etran);
         atran = ptran / (1.0 - etran^2);

         % -----------------  find delta v at point a  ----------------- }
         vinit  = sqrt( mu/rinit );  % sqrt( (2.0 * mu)/rinit - (mu/ainit) );
         vtrana = sqrt( (2.0 * mu)/rinit - (mu/atran) );
         deltava= abs( vtrana - vinit );

         % -----------------  find delta v at point b  ----------------- }
         vfinal  = sqrt( (2.0 * mu)/rfinal - (mu/afinal) );
         vtranb  = sqrt( (2.0 * mu)/rfinal - (mu/atran) );
         fpatranb= atan( ( etran * sin(nutran) ) / ( 1.0 + etran * cos(nutran) ) );
         fpafinal= atan( ( efinal * sin(nutran) ) / ( 1.0 + efinal * cos(nutran) ) );
         deltavb = sqrt( vtranb * vtranb + vfinal * vfinal ...
                          - 2.0 * vtranb * vfinal * cos( fpatranb - fpafinal ) );

         % ----------------  find transfer time of flight  ------------- }
         if etran < 0.99999  
             sinv= ( sqrt( 1.0 - etran * etran ) * sin(nutran) ) ...
                    / ( 1.0 + etran * cos(nutran) );
             cosv= ( etran + cos(nutran) ) / ( 1.0 + etran * cos(nutran) );
             e   = atan22( sinv,cosv );
             dttu= sqrt( (atran * atran * atran)/mu ) * ...
                        ( e - etran * sin(e) - (eainit - etran * sin(eainit)) );
           else
             if abs( etran-1.0 ) < 0.000001  
                 % parabolic dttu }
               else
                 % hyperbolic dttu }
               end; 
           end
       else
         fprintf(1,'the one tangent burn is not possible for this case ' );
       end;

     constastro;
%      fprintf(1,' atran   %11.7f  %11.7f km %11.7f \n',atran, atran*re,ptran*re );
%      fprintf(1,' etran   %11.7f  \n',etran );
%      fprintf(1,' vinit   %11.7f  %11.7f km/s \n',vinit, vinit*velkmps );
%      fprintf(1,' phi  %11.7f  %11.7f km/s \n',fpatranb*rad, fpafinal*rad );
%      fprintf(1,' E  %11.7f  \n',e*rad );
%      fprintf(1,' vtrana  %11.7f  %11.7f km/s \n',vtrana, vtrana*velkmps );
%      fprintf(1,' vtranb  %11.7f  %11.7f km/s \n',vtranb, vtranb*velkmps );
%      fprintf(1,' vfinal  %11.7f  %11.7f km/s \n',vfinal, vfinal*velkmps );

