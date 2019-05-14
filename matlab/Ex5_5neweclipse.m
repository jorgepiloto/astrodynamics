
% STK Sun J2000 Vector on 20 May 2015 18:00:00.000 UTCG:
RSun = [77523937.090045, 119297006.803516, 51716128.096744];
disp('STK Penumbra Start and Stop Eccentric Anomalies: 170.928 and 307.332 deg')
disp('STK Umbra Start and Stop Eccentric Anomalies:    171.456 and 306.803 deg')

% Orbital Elements:
rp  = 6378.1363;
a = 6878.14;
ecc = 0;
incl = 28.5;
raan   = 0;
argp  = 0;
nu  = 0;
mu  = 3.986004417e5;

[Een, Eex] = ShadowEntryExit( RSun, rp, a, ecc, incl, raan, argp, nu, mu );


