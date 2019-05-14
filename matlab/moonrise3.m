% -----------------------------------------------------------------------------
%
%                           function moonriset
%
%  this function finds the universal time for moonrise and moonset given the
%    day and site location.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  revisions
%                - david vallado                                 25 jan 2011
%
%  inputs          description                    range / units
%    jd          - julian date                    days from 4713 bc
%    latgd       - site latitude (south -)        -65ø to 65ø rad
%    lon         - site longitude (west -)        -2pi to 2pi rad
%
%  outputs       :
%    utmoonrise  - universal time of moonrise     hrs
%    utmoonset   - universal time of moonset      hrs
%    moonphaseang- phase angle of the moon        deg
%    error       - error parameter
%
%  locals        :
%    moonangle   - angle between the moon vector
%                  and a point on the earth       rad
%    jdtemp      - julian date for moonrise/set   days from 4713 bc
%    uttemp      - temporary ut time              days
%    tut1        - julian centuries from the
%                  jan 1, 2000 12 h epoch (ut1)
%    rtasc       - right ascension                rad
%    decl        - declination                    rad
%    meanlonmoon -                                rad
%    meananomaly -                                rad
%    eclplong    - longitude of the ecliptic      rad
%    obliquity   - obliquity of the ecliptic      rad
%    rmoon
%    rmoonrs
%    rv
%    rhosat
%    tsry1
%    l, m, n     - direction cosines
%    eclplat
%    moongha, moonghan
%    dgha, dghan
%    lhan
%    lst
%    deltaut, deltautn
%    t, tn
%    hzparal
%    loneclsun
%    loneclmoon
%    ttdb
%    gst         - for 0 h utc of each day        rad
%    lha         - local hour angle               rad
%    year        - year                           1900 .. 2100
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - hour                           0 .. 23
%    min         - minute                         0 .. 59
%    sec         - second                         0.0  .. 59.999
%    opt         - idx to for rise and set calc    1,2
%
%  coupling      :
%    invjulian- finds the year day mon hr min sec from the julian date
%    julianday   - finds the julian date given year, mon day, hr, min, sec
%
%  references    :
%    vallado       2007, 292, Alg 32, Ex 5-4
%
% [utmoonrise,utmoonset,moonphaseang,error] = moonrise( jd,latgd,lon )
% -----------------------------------------------------------------------------

function [utmoonrise,utmoonset,moonphaseang,error] = moonrise3( jd,latgd,lon, show )

    twopi    =  2.0*pi;
    deg2rad  =  pi/180.0;

    % ------------------------  implementation   ------------------
    error    = 'ok';

    % -------------- for once for moonrise (1), ) set (2) ---------
    % -------------- make sure lon is within +- 180 deg -----------
    if ( lon > pi )
        lon = lon - 2.0 * pi;
    end
    if ( lon < -pi )
        lon = lon + twopi;
    end

 %   try1 = 1; % try another approach on the current option
    % opt = 1-2 for rise/set
    for opt = 1 : 2
   		tolerance = 3.5e-4; % Fraction of day: ~30 sec
		GHA = 0.0;
        LHAn = 0.0;
        jdtemp = jd; % assign incoming
		deltaUT = 0.001;
		deltaJD = 0.0;
		loopCount = 1;
        while ((abs(deltaUT) > tolerance) && (loopCount < 10))
            ttdb = ( jdtemp - 2451545.0 ) / 36525.0;
            eclplong = 218.32 + 481267.8813 * ttdb  ...
                + 6.29 * sin( (134.9 + 477198.85 * ttdb) * deg2rad )  ...
                - 1.27 * sin( (259.2 - 413335.38 * ttdb) * deg2rad )  ...
                + 0.66 * sin( (235.7 + 890534.23 * ttdb) * deg2rad )  ...
                + 0.21 * sin( (269.9 + 954397.70 * ttdb) * deg2rad )  ...
                - 0.19 * sin( (357.5 + 35999.05 * ttdb) * deg2rad )   ...
                - 0.11 * sin( (186.6 + 966404.05 * ttdb) * deg2rad );
            eclplat = 5.13 * sin( ( 93.3 + 483202.03 * ttdb) * deg2rad )  ...
                + 0.28 * sin( (228.2 + 960400.87 * ttdb) * deg2rad ) ...
                - 0.28 * sin( (318.3 +  6003.18 * ttdb) * deg2rad )  ...
                - 0.17 * sin( (217.6 - 407332.20 * ttdb) * deg2rad );
            eclplong = rem( eclplong * deg2rad, twopi );
            eclplat  = rem( eclplat * deg2rad, twopi );
            if show == 'y'
                fprintf('%2i %2i%8.5f  %11.7f %11.7f ',opt, loopCount,ttdb, eclplong/deg2rad, eclplat/deg2rad);
            end;
            obliquity = 23.439291 - 0.0130042 * ttdb;
            obliquity = obliquity * deg2rad;
            
            % ------- find the geocentric direction cosines -------
            l = cos( eclplat ) * cos( eclplong );
            m = cos(obliquity) * cos(eclplat) * sin(eclplong) ...
                - sin(obliquity) * sin(eclplat);
            n = sin(obliquity) * cos(eclplat) *  sin(eclplong) ...
                + cos(obliquity) * sin(eclplat);

            rtasc = atan2( m,l );
            % - check that rtasc is in the same quadrant as eclplong
                    if ( eclplong < 0.0  )
                        eclplong = eclplong + twopi;
                    end
            if ( abs( eclplong - rtasc )  >  pi*0.5  )
                rtasc = rtasc + 0.5 * pi * fix( 0.5 + (eclplong - rtasc) / (0.5 *pi) );
            end
            decl = asin( n );

            [lst, gmst] = lstime(lon,jdtemp);

            GHAn = gmst - rtasc;
			LHA = GHAn + lon;
			if (loopCount == 1)
                deltaGHA = 347.81  * deg2rad;
			else
				deltaGHA = ((GHAn - GHA) / deltaUT);
            end

            if ( deltaGHA < 0.0  )
                deltaGHA = deltaGHA + twopi / abs(deltaUT);
            end

            
			cosLHAn = (0.00233 - sin(latgd) * sin(decl)) / (cos(latgd) * cos(decl));
            if show == 'y'
                fprintf(' coslha %8.5f lmn %11.7f %11.7f %11.7f \n',cosLHAn, l, m, n);
            end;

            if (abs(cosLHAn) > 1.0) 
				% No event on this day; advance to the next
				deltaUT = 1;
            if show == 'y'
				fprintf('nothing Advancing one day \n');
            end    
				fprintf('a');            
            else
				LHAn = acos(cosLHAn);
				if (opt == 1)   % only for rise calcs
					LHAn = twopi - LHAn;
                end
				%if (debugging) out.printf("LHAn - LHA = %f\n", LHAn - LHA);
				deltaUT = (LHAn - LHA) / deltaGHA;
				if (deltaUT < -0.5)
					deltaUT = deltaUT + twopi / deltaGHA;
				else 
				    if (deltaUT > 0.5)
					    deltaUT = deltaUT - twopi / deltaGHA;
                    end
                end
				if (deltaJD + deltaUT < 0.0) 
					deltaUT = deltaUT + 1;
	     			fprintf('A'); 
                     if show == 'y'
			         	fprintf('Advancing one day \n');
                     end
                end
				GHA = GHAn;
            end
           
			jdtemp = jdtemp + deltaUT;
			deltaJD = deltaJD + deltaUT;
			loopCount = loopCount + 1;
            if show == 'y'
			     fprintf('%2i %2i %11.7f hrs %8.5f ',opt, loopCount, deltaUT*24, deltaJD*24);
                 fprintf('rtasc %11.7f  decl %8.5f gmst %8.5f GHAn %8.5f dGHA %8.5f  LHAn %8.5f jdtemp %8.5f\n',rtasc/deg2rad, decl/deg2rad, gmst/deg2rad, GHAn/deg2rad, deltaGHA/deg2rad, LHAn/deg2rad, jdtemp);  
            end
    end  % while loop
    
        if opt == 1
            utmoonrise = deltaJD*24;    
        end    
        if opt == 2
            utmoonset = deltaJD*24;    
        end    
     
     end % for opt
   
  if show == 'y'
    fprintf('rise %11.7f  set %11.7f  \n',utmoonrise, utmoonset);     
  end  
    
  %  pause;

    % ------------- determine phase angle of the moon -------------
    meanlong = 280.4606184 + 36000.77005361 * ttdb;
    meanlong = rem( meanlong,360.0 );

    meananomaly = 357.5277233 + 35999.05034 * ttdb;
    meananomaly = rem( meananomaly*deg2rad,twopi );
    if ( meananomaly < 0.0 )
        meananomaly = twopi + meananomaly;
    end

    loneclsun = meanlong + 1.914666471 * sin(meananomaly)  ...
        + 0.019994643 * sin(2.0 * meananomaly);

    loneclmoon = 218.32 + 481267.8813 * ttdb + 6.29 * sin( (134.9 + 477198.85 * ttdb)  ...
        * deg2rad ) - 1.27 * sin( (259.2 -413335.38 * ttdb) * deg2rad ) ...
        + 0.66 * sin( (235.7 + 890534.23 * ttdb) * deg2rad )  ...
        + 0.21 * sin( (269.9 + 954397.70 * ttdb) * deg2rad )  ...
        - 0.19 * sin( (357.5 + 35999.05 * ttdb) * deg2rad )  ...
        - 0.11 * sin( (186.6 + 966404.05 * ttdb) * deg2rad );
    loneclmoon = rem( loneclmoon, 360.0  );

    moonphaseang = loneclmoon - loneclsun;

    if ( moonphaseang < 0.0  )
        moonphaseang = 360.0  + moonphaseang;
    end

