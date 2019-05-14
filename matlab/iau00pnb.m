%
% ----------------------------------------------------------------------------
%
%                           function iau00pnb
%
%  this function calulates the transformation matrix that accounts for the
%    effects of precession-nutation in the iau2000b theory.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  revisions
%    vallado     - consolidate with iau 2000                     14 feb 2005
%
%  inputs          description                    range / units
%    ttt         - julian centuries of tt
%
%  outputs       :
%    nut         - transformation matrix for ire-gcrf
%    deltapsi    - change in longitude            rad
%    l           - delaunay element               rad
%    ll          - delaunay element               rad
%    f           - delaunay element               rad
%    d           - delaunay element               rad
%    omega       - delaunay element               rad
%    many others
%
%  locals        :
%    x           - coordinate                     rad
%    y           - coordinate                     rad
%    s           - coordinate                     rad
%    axs0        - real coefficients for x        rad
%    a0xi        - integer coefficients for x
%    ays0        - real coefficients for y        rad
%    a0yi        - integer coefficients for y
%    ass0        - real coefficients for s        rad
%    a0si        - integer coefficients for s
%    apn         - real coefficients for nutation rad
%    apni        - integer coefficients for nutation
%    appl        - real coefficients for planetary nutation rad
%    appli       - integer coefficients for planetary nutation
%    ttt2,ttt3,  - powers of ttt
%    deltaeps    - change in obliquity            rad
%
%  coupling      :
%    iau00in     - initialize the arrays
%    fundarg     - find the fundamental arguments
%    precess     - find the precession coefficients
%
%  references    :
%    vallado       2004, 212-214
%
% [ deltapsi, pnb, nut, l, l1, f, d, omega, ...
%   lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
% ] = iau00pnb (ttt);
% ----------------------------------------------------------------------------

function [ deltapsi, pnb, prec, nut, l, l1, f, d, omega, ...
           lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
         ] = iau00pnb (ttt);

        sethelp;

        % " to rad
        convrt  = pi / (180.0*3600.0);
        deg2rad = pi / 180.0;

        ttt2 = ttt  * ttt;
        ttt3 = ttt2 * ttt;
        ttt4 = ttt2 * ttt2;
        ttt5 = ttt3 * ttt2;

        % obtain data for calculations form the 2000b theory
        opt = '02';  % a-all, r-reduced, e-1980 theory
        [ l, l1, f, d, omega, ...
          lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
        ] = fundarg( ttt, opt );

        % ---- obtain data coefficients
        [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti] = iau00in;
%        [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, ape, apei, agst, agsti] = iau00in;

        pnsum = 0.0;
        ensum = 0.0;
        for i = 77 : -1 : 1
            tempval = apni(i,1)*l + apni(i,2)*l1 + apni(i,3)*f + apni(i,4)*d + apni(i,5)*omega;
            pnsum = pnsum + (apn(i,1) + apn(i,2)*ttt) * sin(tempval) ...
                          + (apn(i,5) + apn(i,6)*ttt) * cos(tempval);
            ensum = ensum + (apn(i,3) + apn(i,4)*ttt) * cos(tempval) ...
                          + (apn(i,7) + apn(i,8)*ttt) * sin(tempval);
%             pnsum = pnsum + (apn(i,1) + apn(i,2)*ttt) * sin(tempval) ...
%                           + (apn(i,5) ) * cos(tempval);
%             tempval = apei(i,1)*l + apei(i,2)*l1 + apei(i,3)*f + apei(i,4)*d + apei(i,5)*omega;
%             ensum = ensum + (ape(i,3) + ape(i,4)*ttt) * cos(tempval) ...
%                           + (ape(i,7) ) * sin(tempval);
          end;
        % iau2006 approach - does not seem to be correct
        % looks like they still use the iau2000a method and adjust
%        pnsum = 0.0;
%        % data file is not already reveresed
%        for i = 77 : -1 : 1
%            tempval = apni(i,1)*l + apni(i,2)*l1 + apni(i,3)*f + apni(i,4)*d + apni(i,5)*omega;
%            if i > 1320 
%                pnsum = pnsum + (apn(i,1) * sin(tempval) + apn(i,2) * cos(tempval)) * ttt;
%            else
%                pnsum = pnsum + apn(i,1) * sin(tempval) + apn(i,2) * cos(tempval);
%            end;    
%          end;
% 
%        ensum = 0.0;
%        % data file is already reveresed
%        for i = 77 : -1 : 1 
%            tempval = apei(i,1)*l + apei(i,2)*l1 + apei(i,3)*f + apei(i,4)*d + apei(i,5)*omega;
%            if i > 1037 
%                ensum = ensum + (ape(i,1) * cos(tempval) + ape(i,2) * sin(tempval)) * ttt;
%            else
%                ensum = ensum + ape(i,1) * cos(tempval) + ape(i,2) * sin(tempval);
%            end;    
%          end;
%          %  add planetary and luni-solar components.
%        deltapsi = pnsum;  % rad
%        deltaeps = ensum; 

       
       
        % ------ form the planetary arguments
        pplnsum = -0.000135 * convrt;  % " to rad
        eplnsum =  0.000388 * convrt;

        %  add planetary and luni-solar components.
        deltapsi = pnsum + pplnsum;
        deltaeps = ensum + eplnsum;

        [prec,psia,wa,ea,xa] = precess ( ttt, '10' );

        oblo = 84381.406 * convrt; % " to rad
        % or 84381.406????

        % ----------------- find nutation matrix ----------------------
        % mean to true
        a1  = rot1mat(ea + deltaeps);
        a2  = rot3mat(deltapsi);
        a3  = rot1mat(-ea);

        % j2000 to date (precession)
        a4  = rot3mat(-xa);
        a5  = rot1mat(wa);
        a6  = rot3mat(psia);
        a7  = rot1mat(-oblo);

        % icrs to j2000
        a8  = rot1mat(-0.0068192*convrt);
        a9  = rot2mat(0.0417750*sin(oblo)*convrt);
%      a9  = rot2mat(0.0166170*convrt);
        a10 = rot3mat(0.0146*convrt);

        pnb = a10*a9*a8*a7*a6*a5*a4*a3*a2*a1;

prec = a10*a9*a8*a7*a6*a5*a4;

        nut = a3*a2*a1;

        if iauhelp =='y'
            fprintf(1,'p e %11.7f  %11.7f  \n',pnsum*180/pi,ensum*180/pi );
            fprintf(1,'dpsi %11.7f deps %11.7f  \n',deltapsi*180/pi,deltaeps*180/pi );
            fprintf(1,'psia %11.7f wa %11.7f ea %11.7f xa %11.7f  \n',psia*180/pi,wa*180/pi,ea*180/pi,xa*180/pi );
          end;


        % -------------- these are extra not needed for pnb
        if (iaupnhelp == 'y')
            p = psia + ( deltapsi*sin(ea)*cos(xa) - deltaeps*sin(xa) ) / sin(wa);  % rad
            w = wa + deltapsi*sin(ea)*sin(xa) + deltaeps*cos(xa);

            xbar = sin(w)*sin(p);  % rad
            ybar = -sin(oblo)*cos(w) + cos(oblo)*sin(w)*cos(p);

            x = xbar + (-0.0166170 + 0.01460*ybar)*convrt;  % rad
%            x = xbar + (-0.0417750 + 0.01460*ybar)*convrt;  % rad
            y = ybar + (-0.0068192 - 0.01460*xbar)*convrt;

            % -------- now find a
            a = 0.5 + 0.125*(x*x + y*y);

            % -------- now find s
            ssum0 = 0.0;
            for i = 33: -1 : 1
                tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega;
                ssum0 = ssum0 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
              end;
            ssum1 = 0.0;
            for j = 3: -1 : 1
                i = 33 + j;
                tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega;
                ssum1 = ssum1 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
              end;
            ssum2 = 0.0;
            for j = 25: -1 : 1
                i = 33 + 3 + j;
                tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega;
                ssum2 = ssum2 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
              end;
            ssum3 = 0.0;
            for j = 4: -1 : 1
                i = 33 + 3 + 25 + j;
                tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega;
                ssum3 = ssum3 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
              end;
            ssum4 = 0.0;
            for j = 1: -1 : 1
                i = 33 + 3 + 25 + 4 + j;
                tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega;
                ssum4 = ssum4 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
              end;

            s = 0.000094 + 0.00380835*ttt - 0.00011994*ttt2 ...
                - 0.07257409*ttt3 + 0.00002770*ttt4 + 0.00001561*ttt5; % ...
            s = -x*y*0.5 + s*convrt + ssum0 + ssum1*ttt + ssum2*ttt2 + ssum3*ttt3 + ssum4*ttt4;  % rad

            if iauhelp == 'x'
                fprintf(1,'00pnb  x  %14.12f" y  %14.12f" s %14.12f" a %14.12fdeg \n',x/deg2rad*3600,y/deg2rad*3600,s/deg2rad*3600,a/deg2rad );
%                fprintf(1,'p %11.7f w %11.7f xbar %11.7f ybar %11.7f deg \n',p*180/pi,w*180/pi,xbar*180/pi,ybar*180/pi );
              end;

          end;


