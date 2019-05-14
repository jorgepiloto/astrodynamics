    % ------------------------------------------------------------------------------
    %
    %                           procedure combined
    %
    %  this procedure calculates the delta v's for a hohmann transfer for either
    %    circle to circle, or ellipse to circle.
    %
    %  author        : david vallado                  719-573-2600   5 may  2012
    %
    %  inputs          description                    range / units
    %    rinit       - initial position magnitude     er
    %    rfinal      - final position magnitude       er
    %    einit       - eccentricity of first orbit
    %    efinal      - eccentricity of final orbit
    %    nuinit      - true anomaly of first orbit    0 or pi rad
    %    nufinal     - true anomaly of final orbit    0 or pi rad, opp of nuinit
    %
    %  outputs       :
    %    deltava     - change in velocity at point a  er / tu
    %    deltavb     - change in velocity at point b  er / tu
    %    dttu        - time of flight for the trans   tu
    %
    %  locals        :
    %    sme1        - mech energy of first orbit     er2 / tu
    %    sme2        - mech energy of transfer orbit  er2 / tu
    %    sme3        - mech energy of final orbit     er2 / tu
    %    vinit       - velocity of first orbit at a   er / tu
    %    vtransa     - velocity of trans orbit at a   er / tu
    %    vtransb     - velocity of trans orbit at b   er / tu
    %    vfinal      - velocity of final orbit at b   er / tu
    %    ainit       - semimajor axis of first orbit  er
    %    atrans      - semimajor axis of trans orbit  er
    %    afinal      - semimajor axis of final orbit  er
    %
    %  coupling      :
    %    none.
    %
    %  references    :
    %    vallado       2007, 352-359, ex 6-7
    % [deltai1, deltava, deltavb, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );
    %
    % ------------------------------------------------------------------------------

    function [deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2] = combined( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai );

    show = 'n';
    % --------------------  initialize values   ------------------- }
    mu = 1.0; % canonical

    if show == 'y'
        fprintf(1,'rinit %11.7f %11.7f  rfinal %11.7f  %11.7f \n',rinit, rinit*6378.137 , rfinal, rfinal*6378.137  );
    end
    
    ainit  = (rinit * (1.0 + einit * cos(nuinit))) / (1.0 - einit * einit );
    atran  = (rinit + rfinal) * 0.5;
    afinal = (rfinal * (1.0 + efinal * cos(nufinal))) / (1.0 - efinal * efinal );
    sme1 = -mu / (2.0*ainit);
    sme2 = -mu / (2.0*atran);

    if show == 'y'
        fprintf(1,'a1 %11.7f %11.7f  e2 %11.7f \n',a1, a1*6378.137 , e2 );
        fprintf(1,'a2 %11.7f %11.7f \n',a2, a2*6378.137  );
    end

    % -----------------  find delta v at point a  ----------------- }
    vinit  = sqrt( 2.0*( mu/rinit + sme1 ) );
    vtransa= sqrt( 2.0*( mu/rinit + sme2 ) );
    %     fpa2a= atan( ( e2*sin(nu2a) ) / ( 1.0 + e2*cos(nu2a) ) );
    %     fpa1 = atan( ( einit*sin(nuinit) ) / ( 1.0 + einit*cos(nuinit) ) );
    %     deltava= sqrt( vtransa*vtransa + vinit*vinit - 2.0*vtransa*vinit* ...
    %                     ( sin(fpa2a)*sin(fpa1)+cos(fpa2a)*cos(fpa1)*cos(deltai)));

    % -----------------  find delta v at point b  ----------------- }
    vfinal = sqrt( mu/rfinal );  % assumes circular
    vtransb= sqrt( 2.0*( mu/rfinal + sme2 ) );
    %     fpa2b= atan( ( e2*sin(nu2b) ) / ( 1.0 + e2*cos(nu2b) ) );
    %     fpa3 = atan( ( efinal*sin(nufinal) ) / ( 1.0 + efinal*cos(nufinal) ) );

    if show == 'y'
        vkmps = 7.905365719014;
        fprintf(1,'vinit %11.7f %11.7f  vfinal %11.7f  %11.7f \n',vinit, vinit*vkmps, vfinal, vfinal*vkmps );
        fprintf(1,'vtransa %11.7f %11.7f  vtransb %11.7f  %11.7f \n',vtransa, vtransa*vkmps, vtransb, vtransb*vkmps );
    end

    %     deltavb= sqrt( vtransb*vtransb + vfinal*vfinal - 2.0*vtransb*vfinal* ...
    %( sin(fpa2b)*sin(fpa3)+cos(fpa2b)*cos(fpa3)*cos(deltai)));

    % -------------- find proportions of inclination change ---------------
    % ----------------- this is the approximate approach ------------------
    ratio = rfinal/rinit;
    s = 1.0/deltai * atan(sin(deltai)/(ratio^1.5 + cos(deltai) ) );
    if show == 'y'
        fprintf(1,' s %11.7f \n', s );
    end
    deltai1 = s*deltai;
    deltai2 = (1.0-s)*deltai;

    deltava= sqrt( vinit^2  + vtransa^2 - 2.0*vinit*vtransa*cos(deltai1) );
    deltavb= sqrt( vfinal^2 + vtransb^2 - 2.0*vfinal*vtransb*cos(deltai2) );

    dttu= pi * sqrt( (atran * atran * atran)/mu );
    
    % ----------------- figure orientation of the firings -----------------
    gam1 = acos( -(vinit^2+deltava^2-vtransa^2 ) / (2.0*vinit*deltava) );
    gam2 = acos( -(vtransb^2+deltavb^2-vfinal^2 ) / (2.0*vtransb*deltavb) );
    
    



