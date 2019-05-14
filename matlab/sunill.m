% ------------------------------------------------------------------------------
%
%                           function sunill
%
%  this function calculates the illumination due to the sun.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    jd          - julian date                    days
%    lat         - location latitude              rad
%    lon         - location longitdue             rad
%    sunaz       - sun azimuth                    rad
%    sunel       - sun elevation                  rad
%
%  outputs       :
%    sunillum    - sun illumination
%
%  locals        :
%                -
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2001, 295-297, eq 5-9
%
% [sunillum] = sunill   ( jd, lat, lon,sunaz,sunel );
% ------------------------------------------------------------------------------

function [sunillum] = sunill   ( jd, lat, lon,sunaz,sunel );

        deg2rad    =     0.01745329251994;

        % -------------------------  implementation   -----------------
        [rsun,srtasc,sdecl] = sun( jd ); % au's needed for sun ill

        [lst,gst] = lstime( lon,jd );

        lha = lst - srtasc;

        sunel  = asin( sin(sdecl)*sin(lat) + ...
                 cos(sdecl)*cos(lat)*cos(lha) );

        sinv= -sin(lha)*cos(sdecl)*cos(lat)/(cos(sunel)*cos(lat));
        cosv= ( sin(sdecl)-sin(sunel)*sin(lat) )/(cos(sunel)*cos(lat));
        sunaz  = atan2( sinv,cosv );

        sunel= sunel/deg2rad;

        if (sunel > -18.01 )
            x= sunel/90.0;

        if (sunel >= 20)
             l0=  3.74;
             l1=  3.97;
             l2= -4.07;
             l3=  1.47;
           elseif ((sunel >= 5.0) & (sunel < 20.0))
                 l0=   3.05;
                 l1=  13.28;
                 l2= -45.98;
                 l3=  64.33;
               elseif ((sunel >= -0.8) & (sunel < 5.0))
                     l0=    2.88;
                     l1=   22.26;
                     l2= -207.64;
                     l3= 1034.30;
                   elseif ((sunel >= -5.0) & (sunel < -0.8))
                         l0=    2.88;
                         l1=   21.81;
                         l2= -258.11;
                         l3= -858.36;
                       elseif ((sunel >= -12.0) & (sunel < -5.0))
                             l0=    2.70;
                             l1=   12.17;
                             l2= -431.69;
                             l3=-1899.83;
                           elseif ((sunel >= -18.0) & (sunel < -12.0))
                                 l0=   13.84;
                                 l1=  262.72;
                                 l2= 1447.42;
                                 l3= 2797.93;
                               else
                                 l0= 0.0;
                                 l1= 0.0;
                                 l2= 0.0;
                                 l3= 0.0;
                               end

         l1= l0 + l1*x + l2*x*x + l3*x*x*x;
         sunillum= 10.0^l1;

         if ((sunillum < -1e+36) | (sunillum > 999.999 ))
             sunillum= 0.0;
           end
       else
         sunillum= 0.0;
       end

