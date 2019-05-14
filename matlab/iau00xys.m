%
% ----------------------------------------------------------------------------
%
%                           function iau00xys
%
%  this function calulates the transformation matrix that accounts for the
%    effects of precession-nutation in the iau2000 theory.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  revisions
%    vallado     - consolidate with iau 2000                     14 feb 2005
%
%  inputs          description                    range / units
%    ttt         - julian centuries of tt
%    ddx         - eop correction for x           rad
%    ddy         - eop correction for y           rad
%
%  outputs       :
%    nut         - transformation matrix for ire-gcrf
%    x           - coordinate of cip              rad
%    y           - coordinate of cip              rad
%    s           - coordinate                     rad
%
%  locals        :
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
%    l           - delaunay element               rad
%    ll          - delaunay element               rad
%    f           - delaunay element               rad
%    d           - delaunay element               rad
%    omega       - delaunay element               rad
%    deltaeps    - change in obliquity            rad
%    many others
%
%  coupling      :
%    iau00in     - initialize the arrays
%    fundarg     - find the fundamental arguments
%
%  references    :
%    vallado       2004, 212-214
%
% [x,y,s,nut] = iau00xys (ttt, ddx, ddy);
% ----------------------------------------------------------------------------

function [x,y,s,nut] = iau00xys (ttt, ddx, ddy);

        sethelp;

        % " to rad
        convrt  = pi / (180.0*3600.0);
        deg2rad = pi / 180.0;

        ttt2 = ttt  * ttt;
        ttt3 = ttt2 * ttt;
        ttt4 = ttt2 * ttt2;
        ttt5 = ttt3 * ttt2;

        [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti] = iau00in;

        opt = '10';  % a-all, r-reduced, e-1980 theory
        [ l, l1, f, d, omega, ...
          lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate, ...
        ] = fundarg( ttt, opt );

        % ---------------- first find x
        % the iers code puts the constants in here, however
        % don't sum constants in here because they're larger than the last few terms
        xsum0 = 0.0;
        for i = 1306: -1 : 1
            tempval = a0xi(i,1)*l + a0xi(i,2)*l1 + a0xi(i,3)*f + a0xi(i,4)*d + a0xi(i,5)*omega + ...
                      a0xi(i,6)*lonmer  + a0xi(i,7)*lonven  + a0xi(i,8)*lonear  + a0xi(i,9)*lonmar + ...
                      a0xi(i,10)*lonjup + a0xi(i,11)*lonsat + a0xi(i,12)*lonurn + a0xi(i,13)*lonnep + a0xi(i,14)*precrate;
            xsum0 = xsum0 + axs0(i,1)*sin(tempval) + axs0(i,2) * cos(tempval);
          end;
        xsum1 = 0.0;
        % note that the index changes here to j. this is because the a0xi etc
        % indicies go from 1 to 1600, but there are 5 groups. the i index counts through each
        % calculation, and j takes care of the individual summations. note that
        % this same process is used for y and s.
        for j = 253: -1 : 1
            i = 1306 + j;
            tempval = a0xi(i,1)*l + a0xi(i,2)*l1 + a0xi(i,3)*f + a0xi(i,4)*d + a0xi(i,5)*omega + ...
                      a0xi(i,6)*lonmer  + a0xi(i,7)*lonven  + a0xi(i,8)*lonear  + a0xi(i,9)*lonmar + ...
                      a0xi(i,10)*lonjup + a0xi(i,11)*lonsat + a0xi(i,12)*lonurn + a0xi(i,13)*lonnep + a0xi(i,14)*precrate;
            xsum1 = xsum1 + axs0(i,1)*sin(tempval) + axs0(i,2)*cos(tempval);
          end;
        xsum2 = 0.0;
        for j = 36: -1 : 1
            i = 1306 + 253 + j;
            tempval = a0xi(i,1)*l + a0xi(i,2)*l1 + a0xi(i,3)*f + a0xi(i,4)*d + a0xi(i,5)*omega + ...
                      a0xi(i,6)*lonmer  + a0xi(i,7)*lonven  + a0xi(i,8)*lonear  + a0xi(i,9)*lonmar + ...
                      a0xi(i,10)*lonjup + a0xi(i,11)*lonsat + a0xi(i,12)*lonurn + a0xi(i,13)*lonnep + a0xi(i,14)*precrate;
            xsum2 = xsum2 + axs0(i,1)*sin(tempval) + axs0(i,2)*cos(tempval);
          end;
        xsum3 = 0.0;
        for j = 4: -1 : 1
            i = 1306 + 253 + 36 + j;
            tempval = a0xi(i,1)*l + a0xi(i,2)*l1 + a0xi(i,3)*f + a0xi(i,4)*d + a0xi(i,5)*omega + ...
                      a0xi(i,6)*lonmer  + a0xi(i,7)*lonven  + a0xi(i,8)*lonear  + a0xi(i,9)*lonmar + ...
                      a0xi(i,10)*lonjup + a0xi(i,11)*lonsat + a0xi(i,12)*lonurn + a0xi(i,13)*lonnep + a0xi(i,14)*precrate;
            xsum3 = xsum3 + axs0(i,1)*sin(tempval) + axs0(i,2)*cos(tempval);
          end;
        xsum4 = 0.0;
        for j = 1: -1 : 1
            i = 1306 + 253 + 36 + 4 + j;
            tempval = a0xi(i,1)*l + a0xi(i,2)*l1 + a0xi(i,3)*f + a0xi(i,4)*d + a0xi(i,5)*omega + ...
                      a0xi(i,6)*lonmer  + a0xi(i,7)*lonven  + a0xi(i,8)*lonear  + a0xi(i,9)*lonmar + ...
                      a0xi(i,10)*lonjup + a0xi(i,11)*lonsat + a0xi(i,12)*lonurn + a0xi(i,13)*lonnep + a0xi(i,14)*precrate;
            xsum4 = xsum4 + axs0(i,1)*sin(tempval) + axs0(i,2)*cos(tempval);
          end;

        x = -0.01661699 + 2004.19174288*ttt - 0.42721905*ttt2 ...
            -0.19862054*ttt3 - 0.00004605*ttt4 + 0.00000598*ttt5; % "
        x = x*convrt + xsum0 + xsum1*ttt + xsum2*ttt2 + xsum3*ttt3 + xsum4*ttt4;  % rad

        if iauhelp == 'y'
            fprintf(1,'xys x %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n',xsum0/deg2rad,xsum1/deg2rad,xsum2/deg2rad,xsum3/deg2rad,xsum4/deg2rad );
          end;

        % ---------------- now find y
        ysum0 = 0.0;
        for i = 962: -1 : 1
            tempval = a0yi(i,1)*l + a0yi(i,2)*l1 + a0yi(i,3)*f + a0yi(i,4)*d + a0yi(i,5)*omega + ...
                      a0yi(i,6)*lonmer  + a0yi(i,7)*lonven  + a0yi(i,8)*lonear  + a0yi(i,9)*lonmar + ...
                      a0yi(i,10)*lonjup + a0yi(i,11)*lonsat + a0yi(i,12)*lonurn + a0yi(i,13)*lonnep + a0yi(i,14)*precrate;
            ysum0 = ysum0 + ays0(i,1)*sin(tempval) + ays0(i,2) * cos(tempval);
          end;

        ysum1 = 0.0;
        for j = 277: -1 : 1
            i = 962 + j;
            tempval = a0yi(i,1)*l + a0yi(i,2)*l1 + a0yi(i,3)*f + a0yi(i,4)*d + a0yi(i,5)*omega + ...
                      a0yi(i,6)*lonmer  + a0yi(i,7)*lonven  + a0yi(i,8)*lonear  + a0yi(i,9)*lonmar + ...
                      a0yi(i,10)*lonjup + a0yi(i,11)*lonsat + a0yi(i,12)*lonurn + a0yi(i,13)*lonnep + a0yi(i,14)*precrate;
            ysum1 = ysum1 + ays0(i,1)*sin(tempval) + ays0(i,2)*cos(tempval);
          end;
        ysum2 = 0.0;
        for j = 30: -1 : 1
            i = 962 + 277 + j;
            tempval = a0yi(i,1)*l + a0yi(i,2)*l1 + a0yi(i,3)*f + a0yi(i,4)*d + a0yi(i,5)*omega + ...
                      a0yi(i,6)*lonmer  + a0yi(i,7)*lonven  + a0yi(i,8)*lonear  + a0yi(i,9)*lonmar + ...
                      a0yi(i,10)*lonjup + a0yi(i,11)*lonsat + a0yi(i,12)*lonurn + a0yi(i,13)*lonnep + a0yi(i,14)*precrate;
            ysum2 = ysum2 + ays0(i,1)*sin(tempval) + ays0(i,2)*cos(tempval);
          end;
        ysum3 = 0.0;
        for j = 5: -1 : 1
            i = 962 + 277 + 30 + j;
            tempval = a0yi(i,1)*l + a0yi(i,2)*l1 + a0yi(i,3)*f + a0yi(i,4)*d + a0yi(i,5)*omega + ...
                      a0yi(i,6)*lonmer  + a0yi(i,7)*lonven  + a0yi(i,8)*lonear  + a0yi(i,9)*lonmar + ...
                      a0yi(i,10)*lonjup + a0yi(i,11)*lonsat + a0yi(i,12)*lonurn + a0yi(i,13)*lonnep + a0yi(i,14)*precrate;
            ysum3 = ysum3 + ays0(i,1)*sin(tempval) + ays0(i,2)*cos(tempval);
          end;
        ysum4 = 0.0;
        for j = 1: -1 : 1
            i = 962 + 277 + 30 + 5 + j;
            tempval = a0yi(i,1)*l + a0yi(i,2)*l1 + a0yi(i,3)*f + a0yi(i,4)*d + a0yi(i,5)*omega + ...
                      a0yi(i,6)*lonmer  + a0yi(i,7)*lonven  + a0yi(i,8)*lonear  + a0yi(i,9)*lonmar + ...
                      a0yi(i,10)*lonjup + a0yi(i,11)*lonsat + a0yi(i,12)*lonurn + a0yi(i,13)*lonnep + a0yi(i,14)*precrate;
            ysum4 = ysum4 + ays0(i,1)*sin(tempval) + ays0(i,2)*cos(tempval);
          end;

        y = -0.00695078 - 0.02538199*ttt - 22.40725099*ttt2 ...
            +0.00184228*ttt3 + 0.00111306*ttt4 + 0.00000099*ttt5;
        y = y*convrt + ysum0 + ysum1*ttt + ysum2*ttt2 + ysum3*ttt3 + ysum4*ttt4;  % rad

        if iauhelp == 'y'
             fprintf(1,'xys y %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n',ysum0/deg2rad,ysum1/deg2rad,ysum2/deg2rad,ysum3/deg2rad,ysum4/deg2rad );
          end;

        % ---------------- now find s
        ssum0 = 0.0;
        for i = 33: -1 : 1
            tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega + ...
                      a0si(i,6)*lonmer  + a0si(i,7)*lonven  + a0si(i,8)*lonear  + a0si(i,9)*lonmar + ...
                      a0si(i,10)*lonjup + a0si(i,11)*lonsat + a0si(i,12)*lonurn + a0si(i,13)*lonnep + a0si(i,14)*precrate;
            ssum0 = ssum0 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
          end;
        ssum1 = 0.0;
        for j = 3: -1 : 1
            i = 33 + j;
            tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega + ...
                      a0si(i,6)*lonmer  + a0si(i,7)*lonven  + a0si(i,8)*lonear  + a0si(i,9)*lonmar + ...
                      a0si(i,10)*lonjup + a0si(i,11)*lonsat + a0si(i,12)*lonurn + a0si(i,13)*lonnep + a0si(i,14)*precrate;
            ssum1 = ssum1 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
          end;
        ssum2 = 0.0;
        for j = 25: -1 : 1
            i = 33 + 3 + j;
            tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega + ...
                      a0si(i,6)*lonmer  + a0si(i,7)*lonven  + a0si(i,8)*lonear  + a0si(i,9)*lonmar + ...
                      a0si(i,10)*lonjup + a0si(i,11)*lonsat + a0si(i,12)*lonurn + a0si(i,13)*lonnep + a0si(i,14)*precrate;
            ssum2 = ssum2 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
          end;
        ssum3 = 0.0;
        for j = 4: -1 : 1
            i = 33 + 3 + 25 + j;
            tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega + ...
                      a0si(i,6)*lonmer  + a0si(i,7)*lonven  + a0si(i,8)*lonear  + a0si(i,9)*lonmar + ...
                      a0si(i,10)*lonjup + a0si(i,11)*lonsat + a0si(i,12)*lonurn + a0si(i,13)*lonnep + a0si(i,14)*precrate;
            ssum3 = ssum3 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
          end;
        ssum4 = 0.0;
        for j = 1: -1 : 1
            i = 33 + 3 + 25 + 4 + j;
            tempval = a0si(i,1)*l + a0si(i,2)*l1 + a0si(i,3)*f + a0si(i,4)*d + a0si(i,5)*omega + ...
                      a0si(i,6)*lonmer  + a0si(i,7)*lonven  + a0si(i,8)*lonear  + a0si(i,9)*lonmar + ...
                      a0si(i,10)*lonjup + a0si(i,11)*lonsat + a0si(i,12)*lonurn + a0si(i,13)*lonnep + a0si(i,14)*precrate;
            ssum4 = ssum4 + ass0(i,1)*sin(tempval) + ass0(i,2)*cos(tempval);
          end;

        s = 0.000094 + 0.00380835*ttt - 0.00011994*ttt2 ...
            - 0.07257409*ttt3 + 0.00002770*ttt4 + 0.00001561*ttt5; % ...
%            + 0.00000171*ttt*sin(omega) + 0.00000357*ttt*cos(2.0*omega) ...
%            + 0.00074353*ttt2*sin(omega) + 0.00005691*ttt2*sin(2.0*(f-d+omega)) ...
%            + 0.00000984*ttt2*sin(2.0*(f+omega)) - 0.00000885*ttt2*sin(2.0*omega);
        s = -x*y*0.5 + s*convrt + ssum0 + ssum1*ttt + ssum2*ttt2 + ssum3*ttt3 + ssum4*ttt4;  % rad

        if iauhelp == 'y'
            fprintf(1,'xys s %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n',ssum0/deg2rad,ssum1/deg2rad,ssum2/deg2rad,ssum3/deg2rad,ssum4/deg2rad );
          end;

        % add corrections if available
        x = x + ddx;
        y = y + ddy;
          
        % ---------------- now find a
        a = 0.5 + 0.125*(x*x + y*y); % units take on whatever x and y are

        if iauhelp == 'x'
%            fprintf(1,'00xys  x  %14.12f y  %14.12f s %14.12f a %14.12f deg \n',x/deg2rad,y/deg2rad,s/deg2rad,a/deg2rad );
            fprintf(1,'00xys  x  %14.12f" y  %14.12f" s %14.12f" a %14.12fdeg \n',x/deg2rad*3600,y/deg2rad*3600,s/deg2rad*3600,a/deg2rad );
          end;

        % ----------------- find nutation matrix ----------------------
        nut1(1,1) = 1.0 - a*x*x;
        nut1(1,2) = -a*x*y;
        nut1(1,3) = x;
        nut1(2,1) = -a*x*y;
        nut1(2,2) = 1.0 - a*y*y;
        nut1(2,3) = y;
        nut1(3,1) = -x;
        nut1(3,2) = -y;
        nut1(3,3) = 1.0 - a*(x*x + y*y);
%nut1

        nut2 = eye(3);
        nut2(1,1) =  cos(s);
        nut2(2,2) =  cos(s);
        nut2(1,2) =  sin(s);
        nut2(2,1) = -sin(s);

        nut = nut1*nut2;

%       the matrix apears to be orthogonal now, so the extra processing is not needed.
%        if (x ~= 0.0) && (y ~= 0.0)
%            e = atan2(y,x);
%          else
%            e = 0.0;
%          end;
%        d = atan( sqrt((x^2 + y^2) / (1.0-x^2-y^2)) );
%        nut1 = rot3mat(-e)*rot2mat(-d)*rot3mat(e+s)


