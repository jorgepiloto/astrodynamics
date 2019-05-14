%
% ----------------------------------------------------------------------------
%
%                           function iau00in
%
%  this function initializes the matricies needed for iau 2000 reduction
%    calculations. the routine uses the files listed as inputs, but they are
%    are not input to the routine as they are static files.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  revisions
%    dav 14 apr 11  update for iau2011 conventions
%
%  inputs          description                    range / units
%    none
%    iau00x.dat  - file for x coefficient
%    iau00y.dat  - file for y coefficient
%    iau00s.dat  - file for s coefficient
%    iau00n.dat  - file for nutation coefficients
%    iau00pl.dat notused - file for planetary nutation coefficients
%    iau00gs.dat - file for gmst coefficients
%
%  outputs       :
%    axs0        - real coefficients for x        rad
%    a0xi        - integer coefficients for x
%    ays0        - real coefficients for y        rad
%    a0yi        - integer coefficients for y
%    ass0        - real coefficients for s        rad
%    a0si        - integer coefficients for s
%    apn         - real coefficients for nutation rad
%    apni        - integer coefficients for nutation
%    ape         - real coefficients for obliquity rad
%    apei        - integer coefficients for obliquity
%    agst        - real coefficients for gst      rad
%    agsti       - integer coefficients for gst
%
%  locals        :
%    convrt      - conversion factor to radians
%    i           - index
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado     2004, pg 205-219, 910-912
%
% [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti] = iau00in;
% -----------------------------------------------------------------------------

function [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti] = iau00in;
%function [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, ape, apei, agst, agsti] = iau00in;

       % ------------------------  implementation   -------------------
       % " to rad
       convrtu= (0.000001*pi) /(180.0*3600.0);  % if micro arcsecond
       convrtm= (0.001*pi) /(180.0*3600.0);     % if milli arcsecond

       % ------------------------------
       %  note that since all these coefficients have only a single
       %  decimal place, one could store them as integres, and then simply
       %  divide by one additional power of ten. it woul dmake memeory
       %  storage much smaller and potentially faster.
       % ------------------------------

       % xys values
       load iau00x.dat;
       axs0 = iau00x(:,2:3);  % reals
       a0xi = iau00x(:,4:17); % integers
       for i=1:size(axs0)
           axs0(i,1)= axs0(i,1) * convrtu;  % rad
           axs0(i,2)= axs0(i,2) * convrtu;  % rad
       end;

       load iau00y.dat;
       ays0 = iau00y(:,2:3);
       a0yi = iau00y(:,4:17);
       for i=1:size(ays0)
           ays0(i,1)= ays0(i,1) * convrtu;
           ays0(i,2)= ays0(i,2) * convrtu;
       end;

       load iau00s.dat;
       ass0 = iau00s(:,2:3);
       a0si = iau00s(:,4:17);
       for i=1:size(ass0)
           ass0(i,1)= ass0(i,1) * convrtu;
           ass0(i,2)= ass0(i,2) * convrtu;
       end;

       
       
       % nutation values old approach iau2003
       load iau03n.dat;
       apni = iau03n(:,1:5);
       apn  = iau03n(:,7:14);
       for i=1:size(apn)
           apn(i,1)= apn(i,1) * convrtm;
           apn(i,2)= apn(i,2) * convrtm;
           apn(i,3)= apn(i,3) * convrtm;
           apn(i,4)= apn(i,4) * convrtm;
           apn(i,5)= apn(i,5) * convrtm;
           apn(i,6)= apn(i,6) * convrtm;
           apn(i,7)= apn(i,7) * convrtm;
           apn(i,8)= apn(i,8) * convrtm;
       end;

       % planetary nutation values
       load iau03pl.dat;
       appli = iau03pl(:,2:15);
       appl  = iau03pl(:,17:21);  % 21 is extra
       for i=1:size(appl)
           appl(i,1)= appl(i,1) * convrtm;
           appl(i,2)= appl(i,2) * convrtm;
           appl(i,3)= appl(i,3) * convrtm;
           appl(i,4)= appl(i,4) * convrtm;
       end;
       
       
       % nutation values planetary now included new iau2006
%       load iau00n.dat;  % luni-solar
%       apn  = iau00n(:,2:3);
%       apni   = iau00n(:,4:17);
%       for i=1:size(apn)
%           apn(i,1)= apn(i,1) * convrtu;
%           apn(i,2)= apn(i,2) * convrtu;
%       end;
% 
%       load iau00e.dat;  % planetary
%       ape  = iau00n(:,2:3);
%       apei   = iau00n(:,4:17);
%       for i=1:size(ape)
%           ape(i,1)= ape(i,1) * convrtu;
%           ape(i,2)= ape(i,2) * convrtu;
%       end;
       
       % gmst values
       % note - these are very similar to the first 34 elements of iau00s.dat,
       % but they are not the same.
       load iau00gs.dat;
       agst  = iau00gs(:,2:3);
       agsti = iau00gs(:,4:17);
       for i=1:size(agst)
           agst(i,1)= agst(i,1) * convrtu;
           agst(i,2)= agst(i,2) * convrtu;
       end;


