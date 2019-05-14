%
% ----------------------------------------------------------------------------
%
%                           function tod2eci
%
%  this function transforms a vector from the true equator true equinox frame
%    of date (tod), to the mean equator mean equinox (j2000) frame.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%    vallado     - consolidate with iau 2000                     14 feb 2005
%
%  inputs          description                    range / units
%    rtod        - position vector of date
%                    true equator, true equinox   km
%    vtod        - velocity vector of date
%                    true equator, true equinox   km/s
%    atod        - acceleration vector of date
%                    true equator, true equinox   km/s2
%    ttt         - julian centuries of tt         centuries
%
%  outputs       :
%    reci        - position vector eci            km
%    veci        - velocity vector eci            km/s
%    aeci        - acceleration vector eci        km/s2
%
%  locals        :
%    deltapsi    - nutation angle                 rad
%    trueeps     - true obliquity of the ecliptic rad
%    meaneps     - mean obliquity of the ecliptic rad
%    omega       -                                rad
%    nut         - matrix for mod - tod
%
%  coupling      :
%   precess      - rotation for precession        mod - eci
%   nutation     - rotation for nutation          tod - mod
%
%  references    :
%    vallado       2001, 219-220, eq 3-68
%
% [reci,veci,aeci] = tod2eci  ( rtod,vtod,atod,ttt,ddpsi,ddeps );
% ----------------------------------------------------------------------------

function [reci,veci,aeci] = tod2eci  ( rtod,vtod,atod,ttt,ddpsi,ddeps );

        [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

        [deltapsi,trueeps,meaneps,omega,nut] = nutation(ttt,ddpsi,ddeps);

        reci = prec*nut*rtod;

        veci = prec*nut*vtod;

        aeci = prec*nut*atod;


