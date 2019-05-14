%
% -----------------------------------------------------------------------------
%
%                           function iau00gst
%
%  this function finds the iau2000 greenwich sidereal time.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  revisions
%
%  inputs          description                    range / units
%    jdut1       - julian date of ut1             days from 4713 bc
%    ttt         - julian centuries of tt
%    deltapsi    - change in longitude            rad
%    l           - delaunay element               rad
%    ll          - delaunay element               rad
%    f           - delaunay element               rad
%    d           - delaunay element               rad
%    omega       - delaunay element               rad
%    many others for planetary values             rad
%
%  outputs       :
%    gst         - greenwich sidereal time        0 to twopi rad
%    st          - transformation matrix
%
%  locals        :
%    temp        - temporary variable for reals   rad
%    tut1d       - days from the jan 1, 2000 12 h epoch (ut1)
%
%  coupling      :
%    iau00in     - initialize the data arrays
%
%  references    :
%    vallado       2004, 216
%
% [gst,st] = iau00gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
%            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);
% -----------------------------------------------------------------------------

function [gst,st] = iau00gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
                    lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);

        sethelp;

        constastro;

        deg2rad = pi/180.0;
        % " to rad
        convrt  = pi / (180.0*3600.0);

        [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti] = iau00in;

        ttt2 = ttt  * ttt;
        ttt3 = ttt2 * ttt;
        ttt4 = ttt2 * ttt2;
        ttt5 = ttt3 * ttt2;

        % mean obliquity of the ecliptic
        epsa = 84381.448 -   46.84024*ttt - 0.00059*ttt2 + 0.001813*ttt3; % "
        epsa = rem(epsa/3600.0 ,360.0  ); % deg
        epsa = epsa * deg2rad; % rad

        %  evaluate the ee complementary terms
        gstsum0 = 0.0;
        for i = 33: -1 : 1
            tempval = agsti(i,1)*l + agsti(i,2)*l1 + agsti(i,3)*f + agsti(i,4)*d + agsti(i,5)*omega + ...
                      agsti(i,6)*lonmer  + agsti(i,7)*lonven  + agsti(i,8)*lonear  + agsti(i,9)*lonmar + ...
                      agsti(i,10)*lonjup + agsti(i,11)*lonsat + agsti(i,12)*lonurn + agsti(i,13)*lonnep + agsti(i,14)*precrate;
            gstsum0 = gstsum0 + agst(i,1)*sin(tempval) + agst(i,2)*cos(tempval); % rad
        end;
        gstsum1 = 0.0;
        for j = 1: -1 : 1
            i = 33 + j;
            tempval = agsti(i,1)*l + agsti(i,2)*l1 + agsti(i,3)*f + agsti(i,4)*d + agsti(i,5)*omega + ...
                      agsti(i,6)*lonmer  + agsti(i,7)*lonven  + agsti(i,8)*lonear  + agsti(i,9)*lonmar + ...
                      agsti(i,10)*lonjup + agsti(i,11)*lonsat + agsti(i,12)*lonurn + agsti(i,13)*lonnep + agsti(i,14)*precrate;
            gstsum1 = gstsum1 + agst(i,1)*ttt*sin(tempval) + agst(i,2)*ttt*cos(tempval);
        end;

        eect2000 = gstsum0 + gstsum1 * ttt;  % rad

        % equation of the equinoxes
        ee2000 = deltapsi * cos(epsa) + eect2000;  % rad

        %  earth rotation angle
        tut1d= jdut1 - 2451545.0;
        era = twopi * ( 0.7790572732640 + 1.00273781191135448 * tut1d );
        era = rem (era,twopi);  % rad

        %  greenwich mean sidereal time, iau 2000.
        gmst2000 = era + (0.014506 + 4612.15739966*ttt + 1.39667721*ttt2 ...
                   - 0.00009344*ttt3 + 0.00001882*ttt4) * convrt; % " to rad

        gst = gmst2000 + ee2000; % rad

        if iauhelp == 'y'
            fprintf(1,'meanobl %11.7f getsum %11.7f %11.7f eect %11.7f  \n',epsa*180/pi,gstsum0*180/pi,gstsum1*180/pi,eect2000*180/pi );
            fprintf(1,'ee2000 %11.7f gmst2000 %11.7f gst %11.7f  \n',ee2000*180/pi,gmst2000*180/pi,gst*180/pi );
        end;

        % transformation matrix
        st(1,1) =  cos(gst);
        st(1,2) = -sin(gst);
        st(1,3) =  0.0;

        st(2,1) =  sin(gst);
        st(2,2) =  cos(gst);
        st(2,3) =  0.0;

        st(3,1) =  0.0;
        st(3,2) =  0.0;
        st(3,3) =  1.0;



        
        
        