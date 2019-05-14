%     -----------------------------------------------------------------
%
%                              Ex5_5.m
%
%  this file demonstrates example 5-5.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2013
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            12 jun 15  david vallado
%                         original
%  changes :
%            12 jun 15  david vallado
%                         original baseline
%
%     *****************************************************************

constmath;

year = 1994;
mon  =  5;
day  =  20;
hr   =  20;
min  =  0;
sec  =  0.000;
dut1 =  0.0000;
dat  = 0;
xp   =  0.0;
yp   =  0.0;
lod  =  0.0;
timezone= 0;
terms = 2;
order = 106;
[ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
         = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

ttdb

% find coes for Jupiter
a = 5.202603191 + 1.913e-7*ttdb;  % in AU
ecc = 0.04849485 + 0.000163244*ttdb - 4.719e-7*ttdb^2;
incl = (1.303270 - 0.0019872*ttdb + 3.318e-5*ttdb^2 + 9.2e-8 * ttdb^3)/rad;
omega = (100.464441 + 0.1766828*ttdb + 0.000903877*ttdb^2 - 7.032e-6*ttdb^3)/rad;
argp1 = (14.331309 + 0.2155525*ttdb + 0.00072252*ttdb^2 - 4.59e-6*ttdb^3)/rad;
lonmean = (34.351484 + 3034.9056746*ttdb - 0.00008501*ttdb^2 + 4.0e-9*ttdb^3)/rad;

argp1*rad
lonmean*rad

mean = lonmean - argp1;
argp = argp1 - omega;

[eccanom, nu] = newtonm(ecc, mean);

au = 149597870.7;  % km
p = a*(1.0 - ecc*ecc) * au;   % in km

musun = 1.32712428e11;  %  km3/s2
% answer in km/s
[r, v] = coe2rvh ( p, ecc, incl, omega, argp, nu,0.0, 0.0, 0.0, musun );
  
% r in au
r = r/au;
% v in au/day
v = (v/au) *86400;

fprintf(1,'r  %11.6f %11.6f %11.6f  AU \n', r );
fprintf(1,'v  %11.6f %11.6f %11.6f  AU/day \n', v );

fprintf(1,'coes %11.4f %11.4f %11.7f %11.5f %11.5f ', p/au,a,ecc,incl*rad,omega*rad );
fprintf(1,'%11.5f %11.5f %11.5f\n', argp*rad,nu*rad,mean*rad );

eps = 23.440021/rad;
[reci] = rot1 ( r, -eps );
[veci] = rot1 ( v, -eps );

fprintf(1,'reci  %11.6f %11.6f %11.6f  AU \n', reci );
fprintf(1,'veci  %11.6f %11.6f %11.6f  AU/day \n', veci );

% now in km and km/s
fprintf(1,'reci  %11.1f %11.1f %11.1f  km \n', reci*au );
fprintf(1,'veci  %11.5f %11.5f %11.5f  km/s \n', veci*au/86400 );

mag(veci*au/86400)

