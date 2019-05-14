% ------------------------------------------------------------------------------
%
%                           function anglesl
%
%  this function solves the problem of orbit determination using three
%    optical sightings and the method of laplace.
%
%  author        : david vallado                  719-573-2600   24 apr 2003
%
%  inputs          description                    range / units
%    re           - radius earth, sun, etc        km
%    mu           - grav param earth, sun etc     km3/s2
%    tu           - time unit                     sec
%    rtasc1       - right ascension #1            rad
%    rtasc2       - right ascension #2            rad
%    rtasc3       - right ascension #3            rad
%    decl1       - declination #1                rad
%    decl2       - declination #2                rad
%    decl3       - declination #3                rad
%    jd1          - julian date of 1st sighting   days from 4713 bc
%    jd2          - julian date of 2nd sighting   days from 4713 bc
%    jd3          - julian date of 3rd sighting   days from 4713 bc
%    rs1          - ijk site position vector #1   km
%    rs2          - ijk site position vector #2   km
%    rs3          - ijk site position vector #3   km
%
%  outputs        :
%    r            - ijk position vector           km
%    v            - ijk velocity vector           km / s
%
%  locals         :
%    l1           - line of sight vector for 1st
%    l2           - line of sight vector for 2nd
%    l3           - line of sight vector for 3rd
%    ldot         - 1st derivative of l2
%    lddot        - 2nd derivative of l2
%    rs2dot       - 1st derivative of rs2 - vel
%    rs2ddot      - 2nd derivative of rs2
%    t12t13       - (t1-t2) * (t1-t3)
%    t21t23       - (t2-t1) * (t2-t3)
%    t31t32       - (t3-t1) * (t3-t2)
%    i            - index
%    d            -
%    d1           -
%    d2           -
%    d3           -
%    d4           -
%    oldr         - previous iteration on r
%    rho          - range from site to satellite at t2
%    rhodot       -
%    dmat         -
%    d1mat        -
%    d2mat        -
%    d3mat        -
%    d4mat        -
%    earthrate    - angular rotation of the earth
%    l2dotrs      - vector l2 dotted with rs
%    temp         - temporary vector
%    temp1        - temporary vector
%    small        - tolerance
%    roots        -
%
%  coupling       :
%    mag          - magnitude of a vector
%    determinant  - evaluate the determinant of a matrix
%    cross        - cross product of two vectors
%    unit         - unit vector
%    factor       - find the roots of a polynomial
%
%  references     :
%    vallado       2001, 413-417
%
% [r2,v2] = anglesl ( decl1,decl2,decl3,rtasc1,rtasc2,rtasc3,jd1,jd2,jd3,rs1,rs2,rs3 );
% ------------------------------------------------------------------------------

function [r2, v2] = anglesl ( decl1,decl2,decl3,rtasc1,rtasc2, ...
        rtasc3,jd1,jd2,jd3,rs1,rs2,rs3, re, mu, tu );

    % -------------------------  implementation   -------------------------
%     omegaearth   =     0.000072921158553;  % earth rad/s
%     omegaearth   =     0.017202791208627;  % sun rad/s
%     omegaearth   =     2.0 * pi/365.24221897;  % au / day
%     earthrate(1)= 0.0;
%     earthrate(2)= 0.0;
%     earthrate(3)= omegaearth;
    %        tuday        =    58.132440906;
    %        mu           =     1.32712428e11;
    % need to switch these for interplanetary

    %        constant;
    small = 0.00000001;


    % ---------- set middle to 0, find deltas to others -----------
    tau1= (jd1-jd2)*tu;
    tau3= (jd3-jd2)*tu;

    % --------------- find line of sight vectors ------------------
    % should be eci...
    l1(1)= cos(decl1)*cos(rtasc1);
    l1(2)= cos(decl1)*sin(rtasc1);
    l1(3)= sin(decl1);
    l2(1)= cos(decl2)*cos(rtasc2);
    l2(2)= cos(decl2)*sin(rtasc2);
    l2(3)= sin(decl2);
    l3(1)= cos(decl3)*cos(rtasc3);
    l3(2)= cos(decl3)*sin(rtasc3);
    l3(3)= sin(decl3);

    % -------------------------------------------------------------
    %       using lagrange interpolation formula to derive an expression
    %       for l(t), substitute t=t2 and differentiate to obtain the
    %       derivatives of l.
    % -------------------------------------------------------------
    t1t13= 1.0 / (tau1*(tau1-tau3));
    t1t3 = 1.0 / (tau1*tau3);
    t31t3= 1.0 / ((tau3-tau1)*tau3);
    for i= 1 : 3
        ldot(i)=      ( -tau3 * t1t13 )*l1(i) + ...
            ( (-tau1-tau3) * t1t3  )*l2(i) + ...
            ( -tau1 * t31t3 )*l3(i);
        lddot(i)= ( 2.0 * t1t13 )*l1(i) + ...
            ( 2.0 * t1t3  )*l2(i) + ...
            ( 2.0 * t31t3 )*l3(i);
    end
ldot
lddot
    ldotmag = mag( ldot );
    lddotmag = mag( lddot );
    % should these unit vectors use a diff name????????//
    ldot = unit( ldot );
    lddot = unit( lddot );
ldot
lddot

    % ------------------- find 2nd derivative of rs ---------------
%     temp = cross( rs1,rs2 );
%     temp1 = cross( rs2,rs3 );
% 
%     %      needs a different test xxxx%%
%     if ( ( mag(temp) > small ) & ( mag( temp1) > small )  )
%         % ------------ all sightings from one site -----------------
%         %          fix this test here
% these are the same????        
        for i= 1 : 3
            rs2dot(i)=      ( -tau3 * t1t13 )*rs1(i) + ...
                ( (-tau1-tau3) * t1t3  )*rs2(i) + ...
                ( -tau1 * t31t3 )*rs3(i);
            rs2ddot(i)= ( 2.0 * t1t13 )*rs1(i) + ...
                ( 2.0 * t1t3  )*rs2(i) + ...
                ( 2.0 * t31t3 )*rs3(i);
        end

    %    [rs2dot] = cross( earthrate,rs2 );
    %    [rs2ddot] = cross( earthrate,rs2dot );
%     else
%         % ---------- each sighting from a different site ----------
%         for i= 1 : 3
%             rs2dot(i)=      ( -tau3 * t1t13 )*rs1(i) + ...
%                 ( (-tau1-tau3) * t1t3  )*rs2(i) + ...
%                 ( -tau1 * t31t3 )*rs3(i)
%             rs2ddot(i)= ( 2.0 * t1t13 )*rs1(i) + ...
%                 ( 2.0 * t1t3  )*rs2(i) + ...
%                 ( 2.0 * t31t3 )*rs3(i)
%         end
%     end
    rs2dot
    rs2ddot

    for i= 1 : 3
        dmat(i,1) =2.0 * l2(i);
        dmat(i,2) =2.0 * ldot(i);
        dmat(i,3) =2.0 * lddot(i);
        dmat
        % ----------------  position determinants -----------------
        dmat1(i,1) =l2(i);
        dmat1(i,2) =ldot(i);
        dmat1(i,3) =rs2ddot(i);
        dmat2(i,1) =l2(i);
        dmat2(i,2) =ldot(i);
        dmat2(i,3) =rs2(i);

        % ------------  velocity determinants ---------------------
        dmat3(i,1) =l2(i);
        dmat3(i,2) =rs2ddot(i);
        dmat3(i,3) =lddot(i);
        dmat4(i,1) =l2(i);
        dmat4(i,2) =rs2(i);
        dmat4(i,3) =lddot(i);
    end

    dmat1
    dmat2
    d = det(dmat);
    d1= det(dmat1);
    d2= det(dmat2);
    d3= det(dmat3);
    d4= det(dmat4);
    %
    % ---------------  iterate to find rho magnitude ----------------
    %     magr= 1.5   % first guess
    %     writeln( 'input initial guess for magr ' )
    %     readln( magr )
    %     i= 1
    %     repeat
    %         oldr= magr
    %         rho= -2.0*d1/d - 2.0*d2/(magr*magr*magr*d)
    %         magr= sqrt( rho*rho + 2.0*rho*l2dotrs + rs2(4)*rs2(4) )
    %         inc(i)
    %         magr= (oldr - magr ) / 2.0             % simple bissection
    %         writeln( fileout,'rho guesses ',i:2,'rho ',rho:14:7,' magr ',magr:14:7,oldr:14:7 )
    % seems to converge, but wrong numbers
    %         inc(i)
    %     until ( abs( oldr-magr ) < small ) | ( i .ge. 30 )

    if ( abs(d) > 0.000001 )
        % --------------- solve eighth order poly -----------------
        l2dotrs= dot( l2,rs2 );
        poly( 1)=  1.0;  % r2^8th variable%%%%%%%%%%%%%%
        poly( 2)=  0.0;
        poly( 3)=  (l2dotrs*4.0*d1/d - 4.0*d1*d1/(d*d) ...
            - mag(rs2)*mag(rs2) );
        poly( 4)=  0.0;
        poly( 5)=  0.0;
        poly( 6)=  mu*(l2dotrs*4.0*d2/d - 8.0*d1*d2/(d*d) );
        poly( 7)=  0.0;
        poly( 8)=  0.0;
        poly( 9)=  -4.0*mu*d2*d2/(d*d);
        rootarr = roots( poly );
poly
rootarr

        % ------------------ find correct (xx) root ----------------
        bigr2= 0.0;
        for j= 1 : 8
            %                if ( abs( roots(j,2) ) < small )
            %                    writeln( 'root ',j,roots(j,1),' + ',roots(j,2),'j' )
            %        temproot= roots(j,1)*roots(j,1)
            %        temproot= temproot*temproot*temproot*temproot +
            %                  poly(3)*temproot*temproot*temproot + poly(6)*roots(j,1)*temproot + poly(9)
            %                    writeln( fileout,'root ',j,roots(j,1),' + ',roots(j,2),'j  value = ',temproot )
            if ( rootarr(j,1) > bigr2 )
                bigr2= rootarr(j,1);
            end
            %                  end
        end
        fprintf(1,'bigr2 %11.7f ',bigr2);
        fprintf(1,'keep this root ? ');
%        input (bigr2);

        rho= -2.0*d1/d - 2.0*mu*d2 / (bigr2*bigr2*bigr2*d);

        % --------- find the middle position vector ---------------
        for k= 1 : 3
            r2(k)= rho*l2(k) + rs2(k);
        end
        magr2 = mag( r2 );
        % ---------------- find rhodot magnitude ------------------
        rhodot= -d3/d - mu*d4/(magr2*magr2*magr2*d);
        %        writeln( fileout,'rho ',rho:14:7 )
        %        writeln( fileout,'rhodot ',rhodot:14:7 )

        % -------------- find middle velocity vector --------------
        for i= 1 : 3
            v2(i)= rhodot*l2(i) + rho*ldot(i) + rs2dot(i);
        end
    else
        fprintf(1,'determinant value was zero %11.7f ',d );
    end

