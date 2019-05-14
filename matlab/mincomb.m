% ------------------------------------------------------------------------------
%
%                           procedure mincomb
%
%  this procedure calculates the delta v's and the change in inclination
%    necessary for the minimum change in velocity when traveling between two
%    non-coplanar orbits.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    rinit       - initial position magnitude     er
%    rfinal      - final position magnitude       er
%    einit       - ecc of first orbit
%    e2          - ecc of trans orbit
%    efinal      - ecc of final orbit
%    nuinit      - true anomaly of first orbit    0 or pi rad
%    nufinal     - true anomaly of final orbit    0 or pi rad
%    iinit       - incl of the first orbit        rad
%    ifinal      - incl of the second orbit       rad
%
%  outputs       :
%    deltai1     - amount of incl chg req at a    rad
%    deltava     - change in velocity at point a  er / tu
%    deltavb     - change in velocity at point b  er / tu
%    dttu        - time of flight for the trans   tu
%    numiter     - number of iterations
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
%    e2          - eccentricity of second orbit
%
%  coupling      :
%    power       - raise a base to a power
%    asin      - arc sine routine
%
%  references    :
%    vallado       2007, 355, alg 42, table 6-3
%function [deltai,deltai1,deltava,deltavb,dttu ] = mincomb(rinit,rfinal,einit,efinal,nuinit,nufinal,iinit,ifinal);
% ----------------------------------------------------------------------------- }

function [deltai,deltai1,deltava,deltavb,dttu ] = mincomb(rinit,rfinal,einit,efinal,nuinit,nufinal,iinit,ifinal);
     rad  = 57.29577951308230;
     % --------------------  initialize values   -------------------- }
     a1  = (rinit*(1.0+einit*cos(nuinit))) / (1.0 - einit*einit );
     a2  = 0.5 * (rinit+rfinal);
     a3  = (rfinal*(1.0+efinal*cos(nufinal))) / (1.0 - efinal*efinal );
     sme1= -1.0 / (2.0*a1);
     sme2= -1.0 / (2.0*a2);
     sme3= -1.0 / (2.0*a3);

     % ----------- find velocities -------- }
     vinit = sqrt( 2.0*( (1.0/rinit) + sme1 ) );
     v1t   = sqrt( 2.0*( (1.0/rinit) + sme2 ) );

     vfinal= sqrt( 2.0*( (1.0/rfinal) + sme3 ) );
     v3t   = sqrt( 2.0*( (1.0/rfinal) + sme2 ) );

     % ----------- find the optimum change of inclination ----------- }
     tdi = ifinal-iinit;

     temp= (1.0/tdi) * arctan( (power(rfinal/rinit,1.5)-cos(tdi)) / sin(tdi) );
     temp= (1.0/tdi) * arctan( sin(tdi) / (power(rfinal/rinit,1.5)+cos(tdi)) );

     deltava= sqrt( v1t*v1t + vinit*vinit - 2.0*v1t*vinit*cos(temp*tdi) );
     deltavb= sqrt( v3t*v3t + vfinal*vfinal - 2.0*v3t*vfinal*cos(tdi*(1.0-temp)) );

     deltai = temp*tdi;
     deltai1= tdi*(1.0-temp);

     % ----------------  find transfer time of flight  -------------- }
     dttu= pi * sqrt( a2*a2*a2 );

     if show = 'y' then
         dvold= abs(v1t-vinit) + sqrt( v3t*v3t + vfinal*vfinal - 2.0*v3t*vfinal*cos(tdi) );
         fprintf(1,'s = ',temp:11:7,' this uses di in rad ' );
         fprintf(1,'rinit ',rinit:14:7,rinit*6378.137:14:7,' rfinal ',rfinal:14:7,rfinal*6378.137:14:7 );
         fprintf(1,'deltai1 ',deltai*rad:13:7,deltai1*rad:13:7 );
         fprintf(1,'deltava ',deltava:13:7,'deltavb ',deltavb:13:7,' er/tu ' );
         fprintf(1,'deltava ',deltava*7.905365998:13:7,'deltavb ',
                                     deltavb*7.905365998:13:7,' km/s ' );
         fprintf(1,1000*(deltava+deltavb)*7.905365998:13:7,' m/s' );
         fprintf(1,'dv old way ',1000*dvold*7.905365998:13:7,' m/s' );
         fprintf(1,'dttu     ',dttu*13.446851158:13:7,' min' );
       end;

     % ----- iterate to find the optimum change of inclination ----- }
     deltainew  = deltai;         % 1st guess, 0.01 to 0.025 seems good }
     deltai1    = 100.0;              % if going to smaller orbit, should be}
     numiter    = 0;                  % 1.0 - 0.025! }

     while abs(deltainew-deltai1) > 0.000001 do
         deltai1= deltainew;
         deltava= sqrt( v1t*v1t + vinit*vinit - 2.0*v1t*vinit* cos(deltai1) );

         deltavb= sqrt( v3t*v3t + vfinal*vfinal - 2.0*v3t*vfinal* cos(tdi-deltai1) );

         deltainew= asin( (deltava*vfinal*v3t*sin(tdi-deltai1)) / (vinit*v1t*deltavb) );
         inc(numiter);
       end;  % while abs() }
     fprintf(1,'iter di ',deltai1*rad:14:6,'ø',numiter:3,
             (deltava+deltavb)*7905.365998:13:7 );


   end;  % procedure mincomb }

