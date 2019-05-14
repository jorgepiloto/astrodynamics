%
%
% dolam - code insertion for testlam program
%

constastro;
muin = mu;  % km, mum if meters
% ---- do the short way multiple revolution cases ---- }
fprintf( fid,'xx  \n' );
fprintf( 1,'xx 501    psinew         dt       x      a          e \n' );
for i= (nrev)*100 : 500
    dt= i*60.0; % sec
    % make target moving...
    if kepmov == 'y'
        [rtgt1,vtgt1,errork] = kepler( rtgto,vtgto,dt,0 );
    else
        errork = '      ok';
        rtgt1 = rtgto;
        vtgt1 = vtgto;
    end;
    
    if i == nrev*100
        % check min energy condition and min time for that
        cosdeltanu= dot(rinto,rtgt1)/(mag(rinto)*mag(rtgt1));
        chord = sqrt( mag(rinto)^2 + mag(rtgt1)^2 - 2.0*mag(rinto)*mag(rtgt1)*cosdeltanu );
        s     = ( mag(rinto) + mag(rtgt1) + chord )*0.5;
        amin  = s/2.0;
        betam  = 2.0* asin( sqrt( (s-chord)/s ) );  % note that 2*amin is s here
        alpham = pi;
        if ( direc == 's' )
            % this one works - gets min for min a on each case
            ttran = sqrt(amin^3/mu) * (2.0*nrev*pi + alpham -sin(alpham) - betam + sin(betam));
            %    tpar  = (s^1.5 - sin(trangle)*(s-chord)^1.5) * sqrt(2)/(3.0*sqrt(mu));
            %    tpar  =  sqrt(s^3/(8.0*mu)) * (pi - betam + sin(betam));
            tpar  = (s^1.5 - (s-chord)^1.5) * sqrt(2)/(3.0*sqrt(mu));
        else
            ttran = sqrt(amin^3/mu) * (2.0*nrev*pi + alpham -sin(alpham) + betam - sin(betam));
            %    tpar  = (s^1.5 + sin(trangle)*(s-chord)^1.5) * sqrt(2)/(3.0*sqrt(mu));
            %    tpar  =  sqrt(s^3/(8.0*mu)) * (pi + betam - sin(betam));
            % this one works - chk when rest is good if olong/short are +-
            tpar  = (s^1.5 + (s-chord)^1.5) * sqrt(2)/(3.0*sqrt(mu));
        end;
        beta  = 2.0* asin( sqrt( (s-chord)/ (s) ) );  % xxxxxxxxxxxxxx should be 2*a here!!!!
        tmin  = sqrt(amin^3/mu) * ((2.0*nrev+1.0)*pi - beta + sin(beta));
        fprintf(1,'dnu %11.7f mins c %11.7f s %11.7f a %11.7f be %11.7f tranmin %11.7f tmin %11.7f  tpar %11.7f\n', ...
            acos(cosdeltanu)*180/pi, chord, s, amin, betam, ttran, tmin, tpar );
    end;
    
    [vtrans1,vtrans2,errorl] = lambertu( rinto,rtgt1,direc,nrev,dt,fid );
    %          [vtrans1,vtrans2,errorl] = lambertb( rinto,rtgt1,direc,overrev, dt );
    %           if strcmp(errorl, '      ok') ==0 && strcmp(errork, '      ok') ==0
    %               [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (rinto,vtrans1, muin); % of trans orbit
    %               dv1 = mag(vinto - vtrans1);
    %               dv2 = mag(vtrans2 - vtgt1);
    %               fprintf( fid,' %11.5f %11.5f %11.5f %11.5f %11.5f \n',a,ecc,dv1,dv2,dv1+dv2 );
    %           else
    %               fprintf( fid,'  0  0 %s \n',errorl );
    %           end;
end;
