%
% ----------------------------------------------------------------------------
%
%                           function mod2eci
%
%  this function transforms a vector from the mean equator mean equinox of
%    date (mod) to the mean equator mean equinox (j2000) frame.
%
%  author        : david vallado                  719-573-2600   25 jun 2002
%
%  revisions
%    vallado     - consolidate with iau 2000                     14 feb 2005
%
%  inputs          description                    range / units
%    rmod        - position vector of date
%                    mean equator, mean equinox   km
%    vmod        - velocity vector of date
%                    mean equator, mean equinox   km/s
%    amod        - acceleration vector of date
%                    mean equator, mean equinox   km/s2
%    ttt         - julian centuries of tt         centuries
%
%  outputs       :
%    reci        - position vector eci            km
%    veci        - velocity vector eci            km/s
%    aeci        - acceleration vector eci        km/s2
%
%  locals        :
%    none.
%
%  coupling      :
%   precess      - rotation for precession        mod - eci
%
%  references    :
%    vallado       2001, 219-220, eq 3-68
%
% [reci,veci,aeci] = mod2eci  ( rmod,vmod,amod,ttt );
% ----------------------------------------------------------------------------

function [reci,veci,aeci] = mod2eci  ( rmod,vmod,amod,ttt );

        [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

        reci = prec*rmod;

        veci = prec*vmod;

        aeci = prec*amod;

