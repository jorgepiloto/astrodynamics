% ----------------------------------------------------------------------------
%
%                           function iau00i2f
%
%  this function trsnforms a vector from the mean equator mean equniox frame
%    (gcrf), to an earth fixed (itrf) frame.  the results take into account
%    the effects of precession, nutation, sidereal time, and polar motion.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  revisions
%
%  inputs          description                    range / units
%    reci        - position vector eci            km
%    veci        - velocity vector eci            km/s
%    aeci        - acceleration vector eci        km/s2
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    option      - which approach to use          a-2000a, b-2000b, c-2000xys
%    ddx         - eop correction for x           rad
%    ddy         - eop correction for y           rad
%
%  outputs       :
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%
%  locals        :
%    pm          - transformation matrix for itrf-pef
%    st          - transformation matrix for pef-ire
%    nut         - transformation matrix for ire-gcrf
%
%  coupling      :
%   iau00pm      - rotation for polar motion      itrf-pef
%   iau00era     - rotation for earth rotation    pef-ire
%   iau00xys     - rotation for prec/nut          ire-gcrf
%
%  references    :
%    vallado       2004, 205-219
%
% [recef,vecef,aecef] = iau00i2f  ( reci,veci,aeci,ttt,jdut1,lod,xp,yp,option, ddx, ddy );
% ----------------------------------------------------------------------------

function [recef,vecef,aecef] = iau00i2f( reci,veci,aeci,ttt,jdut1,lod,xp,yp,option, ddx, ddy )

  %      sethelp;

        % ---- ceo based, iau2000
        if option == 'c'
            [x,y,s,pnb] = iau00xys (ttt, ddx, ddy);
            [st]  = iau00era (jdut1 );
          end;

        % ---- class equinox based, 2000a
        if option == 'a'
            [ deltapsi, pnb, prec, nut, l, l1, f, d, omega, ...
              lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
            ] = iau00pna (ttt);
            [gst,st] = iau00gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
                       lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);
          end;

        % ---- class equinox based, 2000b
        if option == 'b'
            [ deltapsi, pnb, prec, nut, l, l1, f, d, omega, ...
              lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
            ] = iau00pnb (ttt);
            [gst,st] = iau00gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
                       lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);
          end;

        [pm] = polarm(xp,yp,ttt,'01');

        % ---- setup parameters for velocity transformations
        thetasa= 7.29211514670698e-05 * (1.0  - lod/86400.0 );
        omegaearth = [0; 0; thetasa;];

        rpef  = st'*pnb'*reci;
        recef = pm'*rpef;

        vpef  = st'*pnb'*veci - cross( omegaearth,rpef );
        vecef = pm'*vpef;

        temp  = cross(omegaearth,rpef);
        aecef = pm'*(st'*pnb'*aeci - cross(omegaearth,temp) - 2.0*cross(omegaearth,vpef));

%        if iauhelp == 'y'
            rire  = pnb'*reci;
            vire  = pnb'*veci;
            if (option == 'a') || (option == 'b')
                rmod20  = prec'*reci;
                vmod20  = prec'*veci;
                fprintf(1,'eci           IAU-2006 %c   %14.7f %14.7f %14.7f \n',option, reci );
                fprintf(1,'MOD           IAU-2006 %c   %14.7f %14.7f %14.7f',option, rmod20 );
                fprintf(1,' v %14.9f %14.9f %14.9f',vmod20 );
                fprintf(1,' a %14.9f %14.9f %14.9f\n',aeci );
                fprintf(1,'ERS           IAU-2006 %c   %14.7f %14.7f %14.7f',option, rire );
            end;    
            if option == 'c'
                fprintf(1,'CIRS          IAU-2006 CIO %14.7f %14.7f %14.7f',rire );
            end;    
            fprintf(1,' v %14.9f %14.9f %14.9f\n',vire );
            fprintf(1,'TIRS          IAU-2006 %c   %14.7f %14.7f %14.7f',option, rpef );
            fprintf(1,' v %14.9f %14.9f %14.9f\n',vpef );
%          end;


