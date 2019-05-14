%     -----------------------------------------------------------------
%
%                              Ex7_2.m
%
%  this file demonstrates example 7-2.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%             9 oct 07  david vallado
%                         original
%  changes :
%             9 oct 07  david vallado
%                         original baseline
%
%     *****************************************************************

    constmath;
    constastro;

    re = 6378.137;
    mu = 3.986004418e5;
    tu = 86400.0;

%    re = 6378.145;
%    mu = 3.986005e5;

    %    re = 149597870.0;  % km in 1 au
    %    mu = 1.32712428e11;
    %    tu = 86400.0;

    %    re = 1.0;  % 1 au
    %    mu = 1.0;
    %    tu = 1.0 / 58.132440906; % days in one solar tu

    convrt = pi / (180*3600.0);
    
    casenum = 2;  % 2 book
    
%    typerun = 'l'; % laplace
    typerun = 'g'; % gauss
%    typerun = 'd'; % doubler
%    typerun = 'o'; % gooding

    switch casenum
        case 1  % baseline old book case from stk true topocentric values
            filedat =load('sat1access.dat');
%             % gaussian angles only (vallado book, p. 420 new effort)
%             jd1 = jday(2007, 8, 20, 11, 40, 0.000 );
%             jd2 = jday(2007, 8, 20, 11, 50, 0.000 );
%             %   jd3 = jday(2007, 8, 20, 12,  0, 0.000 );
%             jd3 = jday(2007, 8, 20, 12, 20, 0.000 );  % try unequally spaced obs
% 
%             % topocentric values
%             rtasc1 = -0.4172870/rad;   % right ascension - first sighting
%             rtasc2 = 55.0931551/rad;   % right ascension - second sighting
%             rtasc3 = 134.2826693/rad;  % right ascension - third sighting
%             decl1 = 17.4626616/rad;    % declination - first sighting
%             decl2 = 36.5731946/rad;    % declination - second sighting
%             decl3 = 12.0351097/rad;    % declination - third sighting
% 
%             % new toocentric values from STK 
%             jd1 = jday(2007, 8, 20, 11, 40, 15.000 );
%             jd2 = jday(2007, 8, 20, 11, 50, 15.000 );
%             jd3 = jday(2007, 8, 20, 12, 20, 15.000 );  % try unequally spaced obs
%             rtasc1 = 0.693536/rad;   % right ascension - first sighting
%             rtasc2 = 56.539525/rad;   % right ascension - second sighting
%             rtasc3 = 134.566371/rad;  % right ascension - third sighting
%             decl1 = 18.167962/rad;    % declination - first sighting
%             decl2 = 36.667387/rad;    % declination - second sighting
%             decl3 = 11.826861/rad;    % declination - third sighting
% % point directly in middle works too
% %            jd2 = jday(2007, 8, 20, 12, 0, 15.000 );  % try unequally spaced obs
% %            rtasc2 = 99.515225/rad;  % right ascension - third sighting
% %            decl2 = 30.879333/rad;    % declination - third sighting
% %12 0 15  2979.637847263     8163.337016525     7182.254501689    -5.220599017     3.296266531     0.205034309
% %11 50 15 5897.954130507     5791.046114526     6682.733686585    -4.393910234     4.576816355     1.482423676
%         
%             % site position
%             latgd = 40.0/rad;
%             lon   = -110.0/rad;
%             alt   = 2.0;                     % km

            % at 8-20-07 11:50,
            r2ans = [5897.954130507     5791.046114526     6682.733686585];
            v2ans = [  -4.393910234     4.576816355     1.482423676];
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2ans,v2ans, re, mu);

            year  =  2007;
            mon   =   8;
            day   =  20;
            hr    =  11;
            min   =  50;
            sec   =   15.0000;
            dut1  =  -0.1639883;    
            dat   = 33;
            xp    =  0.210428 * convrt;
            yp    =  0.286899 * convrt;
            lod   =  0.0;
            timezone= 0;
            terms = 0;
            ddpsi = 0.0;
            ddeps = 0.0;
        case 2
            filedat =load('sat11access.dat');
%             % gaussian angles only (vallado book, p. 420 new effort)
%             jd1 = jday(2012, 8, 20, 11, 40, 28.000 );
%             %jd2 = jday(2012, 8, 20, 11, 55, 28.000 );
%             jd2 = jday(2012, 8, 20, 11, 50, 28.000 );% try unequally spaced obs
%             jd3 = jday(2012, 8, 20, 12, 20, 28.000 );  
%             % topocentric values
% %           rtasc1 = 316.520503/rad;   % right ascension - first sighting
% %            rtasc2 = 93.538318/rad;   % right ascension - second sighting
% %            rtasc3 = 140.113543/rad;  % right ascension - third sighting
% %            decl1 =  -14.991625/rad;    % declination - first sighting
% %            decl2 =  32.919464/rad;    % declination - second sighting
% %            decl3 = 7.700143/rad;    % declination - third sighting
%             rtasc1 = 0.939913/rad;   % right ascension - first sighting
%             rtasc2 = 56.883096/rad;   % right ascension - second sighting
%             rtasc3 = 134.823380/rad;  % right ascension - third sighting
%             decl1 =  18.667717/rad;    % declination - first sighting
%             decl2 =  36.930230/rad;    % declination - second sighting
%             decl3 = 11.722965/rad;    % declination - third sighting
%              
%             % site position
%             latgd = 40.0/rad;
%             lon   = -110.0/rad;
%             alt   = 2.0;                     % km

            % at 8-20-12 11:48:28.000 center time,
            r2ans = [6356.48603405     5290.53225676     6511.39697857   ];
            v2ans = [  -4.17294816     4.77654968     1.72027091];
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2ans,v2ans, re, mu);
       
            year  =  2012;
            mon   =   8;
            day   =  20;
            hr    =  11;
            min   =  55;
            sec   =  28.0000;
            dut1  =  -0.6096413;
            dat   = 35; % leap second in July 2012
            xp    =  0.137495 * convrt;
            yp    =  0.342416 * convrt;
            lod   =  0.0;
            timezone= 0;
            terms = 0;
            ddpsi = 0.0;
            ddeps = 0.0;
        case 3
            filedat =load('Sat11Ex1.dat');
        case 4
            filedat =load('Sat11Ex2.dat');
        case 5 
            filedat =load('Sat11Ex3.dat');
        case 6
            filedat =load('Sat11Ex4.dat');
        case 7
            filedat =load('Sat11Ex5.dat');
        case 8
            filedat =load('Sat11Ex6.dat');
        case 9
            filedat =load('Sat11Ex7.dat');
    end; % switch

    
    obs1 = 62; %62 
    obs2 = 70; %70  72
    obs3 = 74; %74 102

    % set eop parameters for new cases
    if casenum >= 3
        dut1  =  0; %-0.6096413;
        dat   = 34; % leap second in July 2012
        xp    =  0.0; %137495 * convrt;
        yp    =  0.0; %342416 * convrt;
        lod   =  0.0;
        timezone= 0;
        terms = 0;
        ddpsi = 0.0;
        ddeps = 0.0;

        obs1 = 1; %62 
        obs2 = 2; %70  72
        obs3 = 3; %74 102
    end;
    
     % ------ read all the data in and process
     numobs = size(filedat,1); % just get # of rows

     %load data into x y z arrays
     yeararr = filedat(:,3); 
     monarr  = filedat(:,2); 
     dayarr  = filedat(:,1); 
     hrarr   = filedat(:,4); 
     minarr  = filedat(:,5); 
     secarr  = filedat(:,6); 
     latarr  = filedat(:,7)/rad; 
     lonarr  = filedat(:,8)/rad; 
     altarr  = filedat(:,9); % km
     rtascarr  = filedat(:,10)/rad; % rad
     declarr   = filedat(:,11)/rad; % rad 

     for j = 1:numobs  % 5 iterations for now
        [obsrecarr(j,1).jd,obsrecarr(j,1).jdf] = jday(yeararr(j),monarr(j),dayarr(j),hrarr(j),minarr(j),secarr(j));
        obsrecarr(j,1).latgd = latarr(j);  % assumes the same sensor site
        obsrecarr(j,1).lon = lonarr(j);  
        obsrecarr(j,1).alt = altarr(j);  
        [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
              = convtime ( yeararr(j), monarr(j), dayarr(j), hrarr(j), minarr(j), secarr(j), 0, dut1, dat );
        [obsrecarr(j,1).rs,obsrecarr(j,1).vs] = site ( latarr(j),lonarr(j),altarr(j) );
        obsrecarr(j,1).ttt = ttt;  
        obsrecarr(j,1).jdut1 = jdut1;
        obsrecarr(j,1).jdut1frac = jdut1frac;
        obsrecarr(j,1).xp = xp;  % rad
        obsrecarr(j,1).yp = yp;  
        obsrecarr(j,1).rtasc = rtascarr(j);
        obsrecarr(j,1).decl = declarr(j);
     end
    

    rtasc1 = obsrecarr(obs1,1).rtasc;  
    rtasc2 = obsrecarr(obs2,1).rtasc;  
    rtasc3 = obsrecarr(obs3,1).rtasc;  
    decl1 = obsrecarr(obs1,1).decl;  
    decl2 = obsrecarr(obs2,1).decl;  
    decl3 = obsrecarr(obs3,1).decl;  
    jd1 = obsrecarr(obs1,1).jd+obsrecarr(obs1,1).jdf;  
    jd2 = obsrecarr(obs2,1).jd+obsrecarr(obs2,1).jdf;  
    jd3 = obsrecarr(obs3,1).jd+obsrecarr(obs3,1).jdf;  
    rs1 = obsrecarr(obs1,1).rs;  % ecef
    vs1 = obsrecarr(obs1,1).vs;
    rs2 = obsrecarr(obs2,1).rs;
    vs2 = obsrecarr(obs2,1).vs;
    rs3 = obsrecarr(obs3,1).rs;
    vs3 = obsrecarr(obs3,1).vs;
    rs1
    
    if casenum < 3
        fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
            p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    end;        
    [year,mon,day,hr,min,second] = invjday(obsrecarr(obs1,1).jd, obsrecarr(obs1,1).jdf);

    utc = second;
    ut1 = utc+dut1;
    tai = utc+dat;
    tt  = tai+32.184;
    [jdut1,jdut1frac] = jday(year,mon,day,hr,min,ut1);
    [jdtt,jdttfrac]  = jday(year,mon,day,hr,min,tt);
    ttt   =  (jdtt-2451545.0)/36525.0;
    fprintf(1,'year %5i ',year);
    fprintf(1,'mon %4i ',mon);
    fprintf(1,'day %3i ',day);
    fprintf(1,'hr %3i:%2i:%8.6f\n',hr,min,second );
    fprintf(1,'dut1 %8.6f s',dut1);
    fprintf(1,' dat %3i s',dat);
    fprintf(1,' xp %8.6f "',xp);
    fprintf(1,' yp %8.6f "',yp);
    fprintf(1,' lod %8.6f s\n',lod);

    % -------------- convert each site vector from ecef to eci -----------------
    a = [0;0;0];   % dummy acceleration variable for the ecef2eci routine
    [year,mon,day,hr,min,sec] = invjday(jd1,0.0);
    [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
        = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
    [rsite1,vseci,aeci] = ecef2eci(rs1,vs1,a,ttt,jdut1+jdut1frac,lod,xp,yp,2,ddpsi,ddeps);

    [year,mon,day,hr,min,sec] = invjday(jd2,0.0);
    [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
        = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
    [rsite2,vseci,aeci] = ecef2eci(rs2,vs2,a,ttt,jdut1+jdut1frac,lod,xp,yp,2,ddpsi,ddeps);

    [year,mon,day,hr,min,sec] = invjday(jd3,0.0);
    [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
        = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
    [rsite3,vseci,aeci] = ecef2eci(rs3,vs3,a,ttt,jdut1+jdut1frac,lod,xp,yp,2,ddpsi,ddeps); % eci


    % ---------------------- run the angles-only routine ------------------
    if typerun == 'l'
        [r2,v2] = anglesl( decl1,decl2,decl3,rtasc1,rtasc2, ...
                          rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3, re, mu, tu );
        processtype = 'anglesl';
    end    
    if typerun == 'd'
        [r2,v2] = anglesdr( decl1,decl2,decl3,rtasc1,rtasc2, ...
                          rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3, re, mu, tu );
        processtype = 'anglesdr';
    end
    if typerun == 'g'
       [r2,v2] = anglesg( decl1,decl2,decl3,rtasc1,rtasc2, ...
                         rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3, re, mu, tu );
        processtype = 'anglesg';
    end
    if typerun == 'o'
       [r2,v2] = anglesgood( decl1,decl2,decl3,rtasc1,rtasc2, ...
                           rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3 );
        processtype = 'anglesgood';
    end

    % -------------- write out answer --------------
    fprintf(1,'\n\ninputs: \n\n');
    [latgc,latgd,lon,alt] = ijk2ll ( rs1 ); % need to use ecef one!!
    fprintf(1,'Site obs1 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite1, latgd*rad, lon*rad, alt*1000 );
    [latgc,latgd,lon,alt] = ijk2ll ( rs2 );
    fprintf(1,'Site obs2 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite2, latgd*rad, lon*rad, alt*1000 );
    [latgc,latgd,lon,alt] = ijk2ll ( rs3 );
    fprintf(1,'Site obs3 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite3, latgd*rad, lon*rad, alt*1000 );
    [year,mon,day,hr,min,sec] = invjday ( jd1, 0.0 );
    fprintf(1,'obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc1*rad, decl1*rad );
    [year,mon,day,hr,min,sec] = invjday ( jd2, 0.0 );
    fprintf(1,'obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc2*rad, decl2*rad );
    [year,mon,day,hr,min,sec] = invjday ( jd3, 0.0 );
    fprintf(1,'Obs#3 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc3*rad, decl3*rad );

    fprintf(1,'\nsolution by %s \n\n', processtype);
    fprintf(1,'r2     %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n',r2/re, r2);
    fprintf(1,'r2 ans %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n',r2ans/re, r2ans);

    fprintf(1,'v2     %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n',v2/velkmps, v2);
    fprintf(1,'v2 ans %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n',v2ans/velkmps, v2ans);

    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2,v2, re, mu);
    fprintf(1,'         p km          a km         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n');
    fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f \n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );

    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2ans,v2ans, re, mu);
    fprintf(1,'         p km          a km         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n');
    fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f \n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );

