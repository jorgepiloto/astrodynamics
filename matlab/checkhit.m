% ------------------------------------------------------------------------------
%
%                           function checkhitearth
%
%  this function checks to see if the trajectory hits the earth during the
%    transfer.  It may calculate quicker if done in canonical units.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    r1          - initial position vector of int km
%    v1t         - initial velocity vector of trnskm/s
%
%  outputs       :
%    hitearth    - is earth was impacted          'y' 'n'
%
%  locals        :
%    sme         - specific mechanical energy
%    rp          - radius of perigee              km
%    transa      - semi-or axis of transfer       km
%    transe      - eccentricity of transfer
%    transp      - semi-paramater of transfer     km
%    hbar        - angular momentum vector of
%                  transfer orbit
%
%  coupling      :
%    dot         - dot product of vectors
%    mag         - magnitude of a vector
%    cross       - cross product of vectors
%
%  references    :
%    vallado       2013, 503, alg 60
%
% [hitearth] = checkhitearth ( rint, v1t );
% ------------------------------------------------------------------------------

function [hitearth] = checkhitearth ( r1,v1t )
        constastro;
        
        % --------------------------  implementation   -----------------
        hitearth= 'n';
        transa = 0.0;
        
        % ----------------------  find h n and e  ----------------------
        hbar = cross( r1,v1t );
        magh = mag( hbar );

        if ( magh > 0.00001  )
            % ---------  find a e and semi-latus rectum   ---------
            sme    = mag(v1t)^2 *0.5  - ( mu /mag(r1) );
            transp = magh*magh/mu;
            if ( abs( sme ) > 0.00001  )
                transa= -mu / (2.0 *sme);
                transe= sqrt( (transa - transp)/transa );
                rp= transa*(1.0 - transe);
            else
                rp= transp*0.5;    % parabola
            end

            if ( abs( rp ) < re + 100.0 )
                hitearth= 'y';
            end
            fprintf( 1, 'hitearth? %s  %11.7f \n',hitearth, transa );
        else
            fprintf( 1,'the orbit does not exist\n');
        end

