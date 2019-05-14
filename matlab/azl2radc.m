
% azl2radc
%
% this function finds the rtasc decl values given the az-el
%
%
%
function [rtasc,decl] = azl2radc(az, el, lat, lst);

    rad = 180.0/pi;

    decl = asin( sin(el)*sin(lat) + cos(el)*cos(lat)*cos(az) );

    slha1 = -( sin(az)*cos(el)*cos(lat) ) / ( cos(decl)*cos(lat) );
    clha1 = ( sin(el) - sin(lat)*sin(decl) ) / ( cos(decl)*cos(lat) );

    lha1 = atan2(slha1,clha1);
%    fprintf(1,' lha1 %13.7f \n',lha1*rad);

    % alt approach
    slha2 = -( sin(az)*cos(el) ) / ( cos(decl) );
    clha2 = ( cos(lat)*sin(el) - sin(lat)*cos(el)*cos(az) ) / ( cos(decl) );

    lha2 = atan2(slha2,clha2);
%    fprintf(1,' lha2 %13.7f \n',lha2*rad);

    rtasc = lst - lha1;


