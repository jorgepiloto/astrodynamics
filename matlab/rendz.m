% ------------------------------------------------------------------------------
%
%                           procedure rendz
%
%  this procedure calculates parameters for a hohmann transfer rendezvous.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    rcs1        - radius of circular orbit int   er
%    rcs2        - radius of circular orbit tgt   er
%    einit       - ecc of first orbit
%    efinal      - ecc of final orbit
%    nuinit      - true anomaly of first orbit    0 or pi rad
%    nufinal     - true anomaly of final orbit    0 or pi rad
%    phasei      - initial phase angle (tgt-int)  +(ahead) or -(behind) rad
%    numrevs     - number of revs to wait
%    ktgt        -
%    kint        -
%
%  outputs       :
%    phasef      - final phase angle              rad
%    waittime    - wait before next intercept opp tu
%    deltav      - change in velocity             er/tu
%
%  locals        :
%    dttutrans   - time of flight of trans orbit  tu
%    atrans      - semimajor axis of trans orbit  er
%    angveltgt   - angular velocity of target     rad / tu
%    angvelint   - angular velocity of int        rad / tu
%    leadang     - lead angle                     rad
%
%  coupling      :
%    power       - raise a base to a power
%
%  references    :
%    vallado       2007, 364, alg 44, alg 45, ex 6-8, ex 6-9
%function [ phasef,waittime,deltav] = rendz(rcs1,rcs3,phasei,einit,efinal,nuinit,nufinal,ktgt,kint);
% ----------------------------------------------------------------------------- }

function [ phasef,waittime,deltav] = rendz(rcs1,rcs3,phasei,einit,efinal,nuinit,nufinal,ktgt,kint);
constastro;
     twopi   =  6.28318530717959;
     mu = 1.0;  % canonical; 
     vkmps = 7.905365719014;

     angvelint = sqrt( mu / (rcs1 * rcs1 * rcs1) );
     angveltgt = sqrt( mu / (rcs3 * rcs3 * rcs3) );
     vint      = sqrt( mu/rcs1 );

  fprintf(1,' angvelint %11.7f %11.7f rad/s  \n',angvelint, angvelint /tusec);
  fprintf(1,' angveltgt %11.7f %11.7f rad/s  \n',angveltgt, angveltgt /tusec);
     
     % ---------- check for satellites in the same orbits ----------- }
     if abs( angvelint - angveltgt ) < 0.000001  
         periodtrans= ( ktgt * twopi + phasei ) / angveltgt;
         atrans     = (periodtrans/(twopi * kint)) ^(2.0/3.0);
         rp         = 2.0 * atrans - rcs1;
         if rp < 1.0 
             fprintf(1,' error - the transfer orbit intersects the earth ' );
           end; 
         vtrans  = sqrt( ((2.0 * mu)/rcs1) - (mu/atrans) );
         deltav  = 2.0 * (vtrans-vint);
 
         phasef= phasei;
         waittime= periodtrans;
         leadang= 0.0;

       fprintf(1,'vint %11.7f %11.7f km/s \n',vint, vint*vkmps  );
       fprintf(1,'vtrans %11.7f %11.7f km/s \n',vtrans, vtrans*vkmps  );
       fprintf(1,'deltavinit %11.7f %11.7f km/s \n',(vtrans-vint), (vtrans-vint)*vkmps  );
       fprintf(1,'atrans %11.7f %11.7f km \n',atrans, atrans*6378.137  );
         
     else
         % ---- different orbits  
         atrans    = (rcs1 + rcs3) / 2.0;
         dttutrans = pi * sqrt( atrans * atrans * atrans / mu );

         leadang = angveltgt * dttutrans;
         phasef  = leadang - pi;
         if phasef < 0.0
              phasef = phasef + pi;
         end
         waittime= ( phasef - phasei + 2.0 * pi * ktgt ) / ( angvelint - angveltgt );

         a1  = (rcs1 * (1.0 + einit * cos(nuinit))) / (1.0 - einit * einit );
         a2  = ( rcs1 + rcs3 ) / 2.0;
         a3  = (rcs3 * (1.0 + efinal * cos(nufinal))) / (1.0 - efinal * efinal );
         sme1= -mu / (2.0 * a1);
         sme2= -mu / (2.0 * a2);
         sme3= -mu / (2.0 * a3);
     % -----------------  find delta v at point a  ------------------ }
         vinit = sqrt( 2.0 * ( (mu/rcs1) + sme1 ) );
         vtransa= sqrt( 2.0 * ( (mu/rcs1) + sme2 ) );
         deltava= abs( vtransa - vinit );

     % -----------------  find delta v at point b  ------------------ }
         vfinal = sqrt( 2.0 * ( (mu/rcs3) + sme3 ) );
         vtransb= sqrt( 2.0 * ( (mu/rcs3) + sme2 ) );
         deltavb= abs( vfinal - vtransb );
         deltav= deltava + deltavb;

       fprintf(1,'leadang %11.7f %11.7f  \n',leadang, leadang * rad);
       fprintf(1,'atrans %11.7f %11.7f km \n',atrans, atrans*6378.137  );
       fprintf(1,'vinit %11.7f %11.7f km/s \n',vinit, vinit*vkmps  );
       fprintf(1,'vtransa %11.7f %11.7f km/s \n',vtransa, vtransa*vkmps  );
       fprintf(1,'vfinal %11.7f %11.7f km/s \n',vfinal, vfinal*vkmps  );
       fprintf(1,'vtransb %11.7f %11.7f km/s \n',vtransb, vtransb*vkmps  );
       fprintf(1,'deltava %11.7f %11.7f km/s \n',deltava, deltava*vkmps  );
       fprintf(1,'deltavb %11.7f %11.7f km/s \n',deltavb, deltavb*vkmps  );

     end;



