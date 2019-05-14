%
% ----------------------------------------------------------------------------
%
%                           function eci2mod
%
%  this function transfroms a vector from the mean equator, mean equinox frame
%    (j2000), to the mean equator mean equinox of date (mod).
%
%  author        : david vallado                  719-573-2600   27 jun 2002
%
%  revisions
%    vallado     - consolidate with iau 2000                     14 feb 2005
%
%  inputs          description                    range / units
%    reci      - position vector eci          km
%    veci      - velocity vector eci          km/s
%    aeci      - acceleration vector eci      km/s2
%    ttt         - julian centuries of tt         centuries
%
%  outputs       :
%    rmod        - position vector of date
%                    mean equator, mean equinox   km
%    vmod        - velocity vector of date
%                    mean equator, mean equinox   km/s
%    amod        - acceleration vector of date
%                    mean equator, mean equinox   km/s2
%
%  locals        :
%    none.
%
%  coupling      :
%   precess      - rotation for precession        eci - mod
%
%  references    :
%    vallado       2001, 214-215, eq 3-57
%
% [prec] = precession  ( ttt );
% ----------------------------------------------------------------------------

function [rmod,vmod,amod] = eci2mod  ( reci,veci,aeci,ttt );

        [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

        rmod=prec'*reci;

        vmod=prec'*veci;

        amod=prec'*aeci;


